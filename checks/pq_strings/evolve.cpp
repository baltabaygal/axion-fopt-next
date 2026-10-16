#include "TempLat/session/sessionguard.h"
#include "pq_model.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace TempLat;

static void numeric(const ParameterParser& p) {
    for (const char* name:{"N","lSide","dt","tMax","t0","a0","H0","omegaEoS","lambda",
                          "tau_origin","tau_star","n_tau","chi_star","core_growth_until","mass_onset","mass_ramp_duration","mass_cap","output_every","snapshot_every"}) {
        const auto at=p.getParams().find(name);if(at==p.getParams().end())continue;
        std::istringstream s(at->second);double v;
        if(!(s>>v)||!std::isfinite(v))throw std::invalid_argument(std::string("Invalid numeric ")+name);
        s>>std::ws;if(!s.eof())throw std::invalid_argument("Trailing numeric text");
        const std::string n=name;
        if((n=="N"||n=="output_every"||n=="snapshot_every")&&(v!=std::floor(v)||v<0.||v>10000000.))
            throw std::invalid_argument("Invalid integer");
    }
}
static std::vector<double> readInitial(const std::string& path,int n,double box) {
    static_assert(std::endian::native==std::endian::little);
    std::ifstream in(path,std::ios::binary);char magic[8];std::uint32_t dimension=0;double length=0;
    in.read(magic,8);in.read(reinterpret_cast<char*>(&dimension),4);in.read(reinterpret_cast<char*>(&length),8);
    if(!in||std::string(magic,8)!="PQIN0001"||dimension!=std::uint32_t(n)||length!=box)
        throw std::invalid_argument("Invalid PQIN0001 initial state/grid");
    std::vector<double> values(4*std::size_t(n)*n*n);
    in.read(reinterpret_cast<char*>(values.data()),values.size()*8);
    if(!in||in.peek()!=std::char_traits<char>::eof())throw std::invalid_argument("Initial state byte count");
    for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite initial state");
    return values;
}
static std::vector<double> state(PQModel& m,int n) {
    auto x=m.fldS(0_c).getFullNDHostView(),y=m.fldS(1_c).getFullNDHostView();
    auto px=m.piS(0_c).getFullNDHostView(),py=m.piS(1_c).getFullNDHostView();
    const auto layout=m.fldS(0_c).getLayout();std::vector<double> v;v.reserve(4*std::size_t(n)*n*n);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        device::IdxArray<3> p;layout.putMemoryIndexFromSpatialLocationInto(p,i,j,k);
        for(int d=0;d<3;++d)if(p[d]<layout.getNGhosts()||p[d]>=layout.getNGhosts()+layout.getLocalSizes()[d])
            throw std::runtime_error("State index outside interior");
        for(double f:{x(p[0],p[1],p[2]),y(p[0],p[1],p[2]),px(p[0],p[1],p[2]),py(p[0],p[1],p[2])}) {
            if(!std::isfinite(f))throw std::runtime_error("Nonfinite evolved field");v.push_back(f);
        }
    }
    return v;
}
static void snapshot(PQModel& m,const std::filesystem::path& path,int n,double box,double u) {
    const auto v=state(m,n);std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PQST0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u,m.conformalTime(u),m.aI,m.chi})out.write(reinterpret_cast<const char*>(&d),8);
    out.write(reinterpret_cast<const char*>(v.data()),v.size()*8);
    if(!out)throw std::runtime_error("Snapshot write failed");
}
static void measure(PQModel& m,std::ofstream& out,int n,double box,double u,bool expansion) {
    const auto v=state(m,n);const double a=m.aI,dx=box/n;
    double xmean=0.,ymean=0.,rmean=0.,rmax=0.,kin=0.,grad=0.,rad=0.,tilt=0.,source=0.,lambdaSource=0.;
    auto index=[n](int i,int j,int k){return 4*((std::size_t(i%n)*n+j%n)*n+k%n);};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        const auto p=index(i,j,k);double x=v[p],y=v[p+1],r=std::hypot(x,y),r2=r*r;
        xmean+=x;ymean+=y;rmean+=r;rmax=std::max(rmax,r);
        kin+=(v[p+2]*v[p+2]+v[p+3]*v[p+3])/(2*std::pow(a,6));
        for(auto q:{index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)})
            grad+=((v[q]-x)*(v[q]-x)+(v[q+1]-y)*(v[q+1]-y))/(2*dx*dx*a*a);
        rad+=.25*m.lambda*(r2-1.)*(r2-1.);tilt+=m.chi*(r2-x*r);
        source+=m.chiDerivative*(r2-x*r);
        lambdaSource+=.25*m.lambdaDerivative*(r2-1.)*(r2-1.);
    }
    const double frequency=std::sqrt(12./(dx*dx)+a*a*(m.lambda*(3*rmax*rmax+1)+6*m.chi));
    if(m.dt*frequency>.5)throw std::runtime_error("dt*field-dependent frequency bound > .5");
    const double count=double(n)*n*n;
    kin/=count;grad/=count;rad/=count;tilt/=count;source/=count;lambdaSource/=count;
    // d rho/d tau = -a'/a*(6K+2G) + explicit partial_tau V.
    const double dilution=expansion?-(6*kin+2*grad)/m.conformalTime(u):0.;
    out<<u<<','<<m.conformalTime(u)<<','<<a<<','<<m.chi<<','<<xmean/count<<','<<ymean/count<<','<<rmean/count
       <<','<<kin<<','<<grad<<','<<rad<<','<<tilt<<','<<kin+grad+rad+tilt<<','<<source+lambdaSource<<','<<dilution
       <<','<<rmax<<','<<m.dt*frequency<<','<<m.lambda<<','<<lambdaSource<<','<<source<<'\n';
}
int main(int argc,char** argv){try{
    SessionGuard guard(argc,argv,false);ParameterParser p(argc,argv);numeric(p);RunParameters<double> run(p);
    const int every=p.get<int>("output_every",10)(),frames=p.get<int>("snapshot_every",0)();
    const double origin=p.get<double>("tau_origin")();
    if(run.N<4||run.N>128||run.eType!=PV2||run.t0!=0.||run.a0!=1.||run.withGWs||run.doResolutionPreserving
       ||run.dt<=0.||run.lSide<=0.||run.dt/run.dx>.25||run.tMax<=0.||every<=0||frames<0)
        throw std::invalid_argument("Require N=4..128, PV2, positive grid/time, dt/dx<=.25");
    if(run.expansion&&(!run.fixedBackground||std::abs(run.omegaEoS-1./3.)>1e-12||std::abs(run.H0-1./origin)>1e-12))
        throw std::invalid_argument("Require fixed radiation background H0=1/tau_origin");
    const double stepsFloat=run.tMax/run.dt;
    if(!std::isfinite(stepsFloat)||stepsFloat>10000000.||std::abs(stepsFloat-std::round(stepsFloat))>1e-8)
        throw std::invalid_argument("Require integral bounded steps");
    const std::filesystem::path output(p.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("result_dir exists");
    auto initial=readInitial(p.get<std::string>("initial_state")(),run.N,run.lSide);
    auto toolbox=MemoryToolBox<3>::makeShared(run.N,1);PQModel model(p,run,toolbox,initial);
    if(model.coreGrowthUntil>0. && !run.expansion)
        throw std::invalid_argument("Core-growth preparation requires radiation expansion");
    const double aEnd=run.expansion?1.+run.tMax/origin:1.;model.setExternalTime(run.tMax);
    // A conservative frequency guard includes tilt curvature and largest lattice mode.
    const double frequency=std::sqrt(12./(run.dx*run.dx)+aEnd*aEnd*(2*model.lambda+6*model.chi));
    if(run.dt*frequency>.5)throw std::invalid_argument("dt*estimated maximum frequency > .5");
    model.setExternalTime(0.);std::filesystem::create_directories(output);
    PositionVerlet<double> integrator(model,run);std::ofstream history(output/"history.csv");
    history<<std::setprecision(17)<<"u,tau,a,chi,x_mean,y_mean,r_mean,kinetic,gradient,radial,tilt,total,source,dilution,r_max,dt_frequency_bound,lambda,lambda_source,chi_source\n";
    snapshot(model,output/"initial.bin",run.N,run.lSide,0.);
    const int steps=std::round(stepsFloat);
    for(int step=0;step<=steps;++step){const double u=step*run.dt;
        if(step%every==0||step==steps){integrator.sync(model,u);measure(model,history,run.N,run.lSide,u,run.expansion);}
        if(frames>0&&step%frames==0){integrator.sync(model,u);std::ostringstream name;
            name<<"frame_"<<std::setw(7)<<std::setfill('0')<<step<<".bin";snapshot(model,output/name.str(),run.N,run.lSide,u);}
        if(step<steps)integrator.evolve(model,u,false);
    }
    snapshot(model,output/"final.bin",run.N,run.lSide,run.tMax);
    std::ofstream parameters(output/"parameters.txt");parameters<<p;
    if(!history||!parameters)throw std::runtime_error("Output write failed");
    std::cout<<"Completed "<<steps<<" PQ PV2 steps: "<<output<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
