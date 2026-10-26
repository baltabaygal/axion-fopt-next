#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>
#include <stdexcept>

struct BirthEvent {double u,x,y,z;};
inline std::vector<BirthEvent> readBirths(const std::string& path,double box,double end) {
    std::vector<BirthEvent> events;if(path.empty()||path=="none")return events;
    std::ifstream in(path);if(!in)throw std::invalid_argument("Cannot open candidate file");
    std::string line;double previous=-1;
    while(std::getline(in,line)) {
        if(line.empty()||line[0]=='#')continue;
        if(events.size()>=10000)throw std::invalid_argument("Candidate budget exceeds 10000");
        std::istringstream row(line);BirthEvent e;
        if(!(row>>e.u>>e.x>>e.y>>e.z))throw std::invalid_argument("Malformed birth event");
        row>>std::ws;if(!row.eof())throw std::invalid_argument("Trailing birth data");
        for(double x:{e.u,e.x,e.y,e.z})if(!std::isfinite(x))throw std::invalid_argument("Nonfinite birth");
        if(e.u<=0||e.u>end||e.u<previous||e.x<0||e.x>=box||e.y<0||e.y>=box||e.z<0||e.z>=box)
            throw std::invalid_argument("Birth time/order/position outside contract");
        previous=e.u;events.push_back(e);
    }
    return events;
}
inline double seedShape(double r,double radius,double width) {
    const double start=radius+2*width,stop=radius+4*width;
    if(r>=stop)return 0;
    double taper=1;
    if(r>start){const double z=(r-start)/(stop-start);taper=1-10*z*z*z+15*z*z*z*z-6*z*z*z*z*z;}
    return .5*(1-std::tanh((r-radius)/width))*taper;
}
struct InsertionCheck {bool accepted;double maxCore;std::size_t cells;};
template<class State>
InsertionCheck inspectCore(const State& v,int n,double box,double a,const BirthEvent& e,double radius,double tolerance) {
    double maximum=0;std::size_t count=0;const double dx=box/n;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        double d2=0;int indices[3]={i,j,k};double centre[3]={e.x,e.y,e.z};
        for(int q=0;q<3;++q){double d=std::abs((indices[q]+.5)*dx-centre[q]);d=std::min(d,box-d);d2+=d*d;}
        if(a*a*d2<=radius*radius){++count;maximum=std::max(maximum,std::abs(v[4*((std::size_t(i)*n+j)*n+k)]));}
    }
    return {count>0&&maximum<=tolerance,maximum,count};
}
inline void insertBirth(TempLat::BubbleModel& m,int n,double box,double a,const BirthEvent& e,double radius,double width) {
    using namespace TempLat;
    auto s=m.fldS(0_c).getFullNDHostView(),p=m.piS(0_c).getFullNDHostView();const auto layout=m.fldS(0_c).getLayout();
    const double dx=box/n;double centre[3]={e.x,e.y,e.z};
    for(ptrdiff_t i=0;i<ptrdiff_t(s.extent(0));++i)for(ptrdiff_t j=0;j<ptrdiff_t(s.extent(1));++j)for(ptrdiff_t k=0;k<ptrdiff_t(s.extent(2));++k){
        device::IdxArray<3> xyz;layout.putSpatialLocationFromMemoryIndexInto(xyz,i,j,k);double r2=0;
        for(int d=0;d<3;++d){const auto ix=(xyz[d]%n+n)%n;double distance=std::abs((ix+.5)*dx-centre[d]);distance=std::min(distance,box-distance);r2+=distance*distance;}
        const double b=seedShape(a*std::sqrt(r2),radius,width);
        s(i,j,k)=s(i,j,k)+(1-s(i,j,k))*b;p(i,j,k)=(1-b)*p(i,j,k);
    }
    m.fldS(0_c).getMemoryManager()->pushHostView();m.piS(0_c).getMemoryManager()->pushHostView();
    m.fldS(0_c).setGhostsAreStale();m.piS(0_c).setGhostsAreStale();
}
