#include "TempLat/session/sessionguard.h"
#include "pocket_model.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace TempLat;

static void numeric(const ParameterParser& p) {
    for (const char* name:{"N","lSide","dt","tMax","lambda","delta_v","mass_tv",
                          "output_every","snapshot_every","slice_every"}) {
        const auto at=p.getParams().find(name);if(at==p.getParams().end())continue;
        std::istringstream s(at->second);double v;
        if(!(s>>v)||!std::isfinite(v))throw std::invalid_argument(std::string("Invalid numeric ")+name);
        s>>std::ws;if(!s.eof())throw std::invalid_argument("Trailing numeric text");
        const std::string n=name;
        if((n=="N"||n=="output_every"||n=="snapshot_every"||n=="slice_every")&&(v!=std::floor(v)||v<0.||v>10000000.))
            throw std::invalid_argument("Invalid integer");
    }
}
static std::vector<double> readInitial(const std::string& path,int n,double box) {
    static_assert(std::endian::native==std::endian::little);
    std::ifstream in(path,std::ios::binary);char magic[8];std::uint32_t dimension=0;double length=0;
    in.read(magic,8);in.read(reinterpret_cast<char*>(&dimension),4);in.read(reinterpret_cast<char*>(&length),8);
    if(!in||std::string(magic,8)!="PKIN0001"||dimension!=std::uint32_t(n)||length!=box)
        throw std::invalid_argument("Invalid PKIN0001 initial state/grid");
    std::vector<double> values(4*std::size_t(n)*n*n);
    in.read(reinterpret_cast<char*>(values.data()),values.size()*8);
    if(!in||in.peek()!=std::char_traits<char>::eof())throw std::invalid_argument("Initial state byte count");
    for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite initial state");
    return values;
}
// Cell order (x,y,z) with z fastest; per cell phi, a, pi_phi, pi_a.
static std::vector<double> state(PocketModel& m,int n) {
    auto phi=m.fldS(0_c).getFullNDHostView(),ax=m.fldS(1_c).getFullNDHostView();
    auto pphi=m.piS(0_c).getFullNDHostView(),pax=m.piS(1_c).getFullNDHostView();
    const auto layout=m.fldS(0_c).getLayout();std::vector<double> v;v.reserve(4*std::size_t(n)*n*n);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        device::IdxArray<3> p;layout.putMemoryIndexFromSpatialLocationInto(p,i,j,k);
        for(int d=0;d<3;++d)if(p[d]<layout.getNGhosts()||p[d]>=layout.getNGhosts()+layout.getLocalSizes()[d])
            throw std::runtime_error("State index outside interior");
        for(double f:{phi(p[0],p[1],p[2]),ax(p[0],p[1],p[2]),pphi(p[0],p[1],p[2]),pax(p[0],p[1],p[2])}) {
            if(!std::isfinite(f))throw std::runtime_error("Nonfinite evolved field");v.push_back(f);
        }
    }
    return v;
}
static void snapshot(const std::vector<double>& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKST0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    out.write(reinterpret_cast<const char*>(v.data()),v.size()*8);
    if(!out)throw std::runtime_error("Snapshot write failed");
}
// Plane through the box centre, z index n/2, all four components.
static void slice(const std::vector<double>& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKSL0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        out.write(reinterpret_cast<const char*>(&v[4*((std::size_t(i)*n+j)*n+n/2)]),32);
    if(!out)throw std::runtime_error("Slice write failed");
}
static void measure(PocketModel& m,const std::vector<double>& v,std::ofstream& out,int n,double box,double u,double dt) {
    const double dx=box/n,cell=dx*dx*dx;
    double kd=0.,gd=0.,vd=0.,ka=0.,ga=0.,coupling=0.,fv=0.,axFv=0.,axWall=0.,axTv=0.,phi2Fv=0.,aMax=0.,curvature=0.;
    double c[3]={0.,0.,0.},s[6]={0.,0.,0.,0.,0.,0.};
    auto index=[n](int i,int j,int k){return 4*((std::size_t(i%n)*n+j%n)*n+k%n);};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        const auto p=index(i,j,k);const double phi=v[p],a=v[p+1];
        double gPhi=0.,gA=0.;
        for(auto q:{index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)}){
            gPhi+=(v[q]-phi)*(v[q]-phi);gA+=(v[q+1]-a)*(v[q+1]-a);
        }
        gPhi/=2*dx*dx;gA/=2*dx*dx;
        const double kPhi=.5*v[p+2]*v[p+2],kA=.5*v[p+3]*v[p+3],c2=.5*m.mass2*phi*phi*a*a;
        const double dil=pocketValue<0>(phi,a,m.lambda,m.deltaV,m.mass2)-c2;
        kd+=kPhi;gd+=gPhi;vd+=dil;ka+=kA;ga+=gA;coupling+=c2;
        const double axion=kA+gA+c2;
        if(phi<.5){
            fv+=1.;axFv+=axion;phi2Fv+=phi*phi;
            const double x[3]={(i+.5)*dx,(j+.5)*dx,(k+.5)*dx};
            for(int d=0;d<3;++d)c[d]+=x[d];
            s[0]+=x[0]*x[0];s[1]+=x[1]*x[1];s[2]+=x[2]*x[2];s[3]+=x[0]*x[1];s[4]+=x[0]*x[2];s[5]+=x[1]*x[2];
        } else if(phi<=.99) axWall+=axion; else axTv+=axion;
        aMax=std::max(aMax,std::abs(a));
        const double vPhiPhi=.5*m.lambda*(1.-6.*phi+6.*phi*phi)+m.deltaV*(36.*phi*phi-24.*phi)+m.mass2*a*a;
        curvature=std::max({curvature,vPhiPhi,m.mass2*phi*phi});
    }
    const double frequency=std::sqrt(12./(dx*dx)+curvature);
    if(dt*frequency>.5)throw std::runtime_error("dt*field-dependent frequency bound > .5");
    double centroid[3]={0.,0.,0.},second[6]={0.,0.,0.,0.,0.,0.};
    if(fv>0.){
        for(int d=0;d<3;++d)centroid[d]=c[d]/fv;
        const int pairs[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
        for(int e=0;e<6;++e)second[e]=s[e]/fv-centroid[pairs[e][0]]*centroid[pairs[e][1]];
    }
    const double volume=fv*cell,radius=std::cbrt(3.*volume/(4.*M_PI));
    out<<u<<','<<(kd+gd+vd+ka+ga+coupling)*cell<<','<<kd*cell<<','<<gd*cell<<','<<vd*cell<<','
       <<ka*cell<<','<<ga*cell<<','<<coupling*cell<<','<<volume<<','<<radius<<','
       <<axFv*cell<<','<<axWall*cell<<','<<axTv*cell<<','<<(fv>0.?phi2Fv/fv:0.)<<','<<aMax;
    for(double x:centroid)out<<','<<x;
    for(double x:second)out<<','<<x;
    out<<','<<dt*frequency<<'\n';
}
int main(int argc,char** argv){try{
    SessionGuard guard(argc,argv,false);ParameterParser p(argc,argv);numeric(p);RunParameters<double> run(p);
    const int every=p.get<int>("output_every",10)(),frames=p.get<int>("snapshot_every",0)(),
              slices=p.get<int>("slice_every",0)();
    if(run.N<4||run.N>128||run.eType!=PV2||run.t0!=0.||run.a0!=1.||run.withGWs||run.doResolutionPreserving
       ||run.expansion||run.dt<=0.||run.lSide<=0.||run.dt/run.dx>.25||run.tMax<=0.||every<=0||frames<0||slices<0)
        throw std::invalid_argument("Require N=4..128, PV2, static background, positive grid/time, dt/dx<=.25");
    const double stepsFloat=run.tMax/run.dt;
    if(!std::isfinite(stepsFloat)||stepsFloat>10000000.||std::abs(stepsFloat-std::round(stepsFloat))>1e-8)
        throw std::invalid_argument("Require integral bounded steps");
    const std::filesystem::path output(p.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("result_dir exists");
    auto initial=readInitial(p.get<std::string>("initial_state")(),run.N,run.lSide);
    auto toolbox=MemoryToolBox<3>::makeShared(run.N,1);PocketModel model(p,run,toolbox,initial);
    std::filesystem::create_directories(output);
    PositionVerlet<double> integrator(model,run);std::ofstream history(output/"history.csv");
    history<<std::setprecision(17)<<"u,total,dilaton_kinetic,dilaton_gradient,dilaton_potential,"
           "axion_kinetic,axion_gradient,coupling,fv_volume,fv_radius,axion_fv,axion_wall,axion_tv,"
           "phi2_fv,a_max,cx,cy,cz,sxx,syy,szz,sxy,sxz,syz,dt_frequency_bound\n";
    const int steps=std::round(stepsFloat);
    for(int step=0;step<=steps;++step){const double u=step*run.dt;
        const bool record=step%every==0||step==steps,frame=frames>0&&step%frames==0,plane=slices>0&&step%slices==0;
        if(record||frame||plane||step==0||step==steps){
            integrator.sync(model,u);const auto v=state(model,run.N);
            if(record)measure(model,v,history,run.N,run.lSide,u,run.dt);
            if(step==0)snapshot(v,output/"initial.bin",run.N,run.lSide,u);
            if(step==steps)snapshot(v,output/"final.bin",run.N,run.lSide,u);
            std::ostringstream name;name<<std::setw(7)<<std::setfill('0')<<step<<".bin";
            if(frame)snapshot(v,output/("frame_"+name.str()),run.N,run.lSide,u);
            if(plane)slice(v,output/("slice_"+name.str()),run.N,run.lSide,u);
        }
        if(step<steps)integrator.evolve(model,u,false);
    }
    std::ofstream parameters(output/"parameters.txt");parameters<<p;
    if(!history||!parameters)throw std::runtime_error("Output write failed");
    std::cout<<"Completed "<<steps<<" pocket PV2 steps: "<<output<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
