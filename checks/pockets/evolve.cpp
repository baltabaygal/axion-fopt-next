#include "TempLat/session/sessionguard.h"
#include "pocket_model.h"
#include "pocket_monitor.h"
#include "tracking_monitor.h"
#include "initial_state.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace TempLat;

static void numeric(const ParameterParser& p) {
    for (const char* name:{"N","lSide","dt","tMax","lambda","delta_v","mass_tv","mass_power","vacuum_stiffening","mass_fv","mass_kappa","axion_fa","far_margin","monitor_radius","tracking_inner","tracking_outer",
                          "output_every","snapshot_every","slice_every"}) {
        const auto at=p.getParams().find(name);if(at==p.getParams().end())continue;
        std::istringstream s(at->second);double v;
        if(!(s>>v)||!std::isfinite(v))throw std::invalid_argument(std::string("Invalid numeric ")+name);
        s>>std::ws;if(!s.eof())throw std::invalid_argument("Trailing numeric text");
        const std::string n=name;
        if((n=="N"||n=="mass_power"||n=="output_every"||n=="snapshot_every"||n=="slice_every")&&(v!=std::floor(v)||v<0.||v>10000000.))
            throw std::invalid_argument("Invalid integer");
    }
}
// Cell order (x,y,z) with z fastest; per cell phi, a, pi_phi, pi_a.
template<class View,class Layout>
class PocketStateView {
    std::array<View,4> fields_;
    Layout layout_;
    int n_;
public:
    PocketStateView(std::array<View,4> fields,Layout layout,int n):fields_(fields),layout_(layout),n_(n){}
    std::size_t size() const {return 4*std::size_t(n_)*n_*n_;}
    double operator[](std::size_t scalar) const {
        const auto cell=scalar/4;
        const int i=cell/(std::size_t(n_)*n_),j=(cell/n_)%n_,k=cell%n_;
        device::IdxArray<3> p;layout_.putMemoryIndexFromSpatialLocationInto(p,i,j,k);
        for(int d=0;d<3;++d)if(p[d]<layout_.getNGhosts()||p[d]>=layout_.getNGhosts()+layout_.getLocalSizes()[d])
            throw std::runtime_error("State index outside interior");
        return fields_[scalar%4](p[0],p[1],p[2]);
    }
};
static auto state(PocketModel& m,int n) {
    auto phi=m.fldS(0_c).getFullNDHostView(),ax=m.fldS(1_c).getFullNDHostView();
    auto pphi=m.piS(0_c).getFullNDHostView(),pax=m.piS(1_c).getFullNDHostView();
    const auto layout=m.fldS(0_c).getLayout();
    PocketStateView v(std::array{phi,ax,pphi,pax},layout,n);
    for(std::size_t scalar=0;scalar<v.size();++scalar)
        if(!std::isfinite(v[scalar]))throw std::runtime_error("Nonfinite evolved field");
    return v;
}
template<class State>
static void snapshot(const State& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKST0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    std::array<double,8192> buffer;
    for(std::size_t start=0;start<v.size();start+=buffer.size()){
        const auto count=std::min(buffer.size(),v.size()-start);
        for(std::size_t j=0;j<count;++j)buffer[j]=v[start+j];
        out.write(reinterpret_cast<const char*>(buffer.data()),count*8);
    }
    if(!out)throw std::runtime_error("Snapshot write failed");
}
// Plane through the box centre, z index n/2, all four components.
template<class State>
static void slice(const State& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKSL0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j){
        const auto p=4*((std::size_t(i)*n+j)*n+n/2);
        const std::array<double,4> cell={v[p],v[p+1],v[p+2],v[p+3]};
        out.write(reinterpret_cast<const char*>(cell.data()),32);
    }
    if(!out)throw std::runtime_error("Slice write failed");
}
template<class State>
static void measure(PocketModel& m,const State& v,std::ofstream& out,int n,double box,double u,double dt,double margin,const PocketMonitor& monitor,const TrackingMonitor& tracking) {
    const double dx=box/n,cell=dx*dx*dx;
    AccurateSum kd,gd,vd,ka,ga,coupling,axFv,axWall,axTv,axFar;
    double fv=0.,phi2Fv=0.,aMax=0.,curvature=0.;
    double c[3]={0.,0.,0.},s[6]={0.,0.,0.,0.,0.,0.};
    auto index=[n](int i,int j,int k){return 4*((std::size_t(i%n)*n+j%n)*n+k%n);};
    // First pass: false-vacuum centroid and volume radius, for the far-field sphere.
    double count=0.,mean[3]={0.,0.,0.};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k)if(v[index(i,j,k)]<.5){
        count+=1.;mean[0]+=(i+.5)*dx;mean[1]+=(j+.5)*dx;mean[2]+=(k+.5)*dx;
    }
    if(count>0.)for(double& x:mean)x/=count;
    const double far=std::cbrt(3.*count*cell/(4.*M_PI))+margin;
    auto wrapped=[box](double d){d=std::abs(d);return std::min(d,box-d);};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        const auto p=index(i,j,k);const double phi=v[p],a=v[p+1];
        double gPhi=0.,gA=0.;
        for(auto q:{index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)}){
            gPhi+=(v[q]-phi)*(v[q]-phi);gA+=(v[q+1]-a)*(v[q+1]-a);
        }
        gPhi/=2*dx*dx;gA/=2*dx*dx;
        const double phi2=phi*phi, massShape=m.massPower==4 ? phi2*phi2 : phi2;
        const double kPhi=.5*v[p+2]*v[p+2],kA=.5*v[p+3]*v[p+3],c2=m.useFormation?pockets::formationPotential(phi,a,m.formation).axion:.5*m.mass2*massShape*a*a;
        const double dil=pocketValue<0>(phi,a,m.lambda,m.deltaV,m.mass2,m.massPower,m.vacuumStiffening,m.formationParameters())-c2;
        kd+=kPhi;gd+=gPhi;vd+=dil;ka+=kA;ga+=gA;coupling+=c2;
        const double axion=kA+gA+c2;
        if(phi<.5){
            fv+=1.;axFv+=axion;phi2Fv+=phi*phi;
            const double x[3]={(i+.5)*dx,(j+.5)*dx,(k+.5)*dx};
            for(int d=0;d<3;++d)c[d]+=x[d];
            s[0]+=x[0]*x[0];s[1]+=x[1]*x[1];s[2]+=x[2]*x[2];s[3]+=x[0]*x[1];s[4]+=x[0]*x[2];s[5]+=x[1]*x[2];
        } else if(phi<=.99) axWall+=axion; else axTv+=axion;
        if(count>0.){
            const double d0=wrapped((i+.5)*dx-mean[0]),d1=wrapped((j+.5)*dx-mean[1]),d2=wrapped((k+.5)*dx-mean[2]);
            if(d0*d0+d1*d1+d2*d2>far*far)axFar+=axion;
        }
        aMax=std::max(aMax,std::abs(a));
        double vPhiPhi=pocketDilatonCurvature(phi,m.lambda,m.deltaV,m.vacuumStiffening)
            +m.mass2*a*a*(m.massPower==4 ? 6.*phi2 : 1.);
        double vAA=m.mass2*massShape;
        double vPhiA=m.mass2*a*(m.massPower==4 ? 4.*phi2*phi : 2.*phi);
        if(m.useFormation) {
            const auto h=pockets::formationPotential(phi,a,m.formation);
            vPhiPhi=h.hphi_phi;vAA=h.ha_a;vPhiA=h.hphi_a;
        }
        const double largest=.5*(vPhiPhi+vAA+std::hypot(vPhiPhi-vAA,2.*vPhiA));
        curvature=std::max(curvature,largest);
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
    std::array<double,21> tracked{};
    if(tracking.enabled())tracked=tracking.measure(v,m.lambda,m.deltaV,m.mass2,m.massPower,m.vacuumStiffening,m.formationParameters());
    out<<u<<','<<(kd+gd+vd+ka+ga+coupling)*cell<<','<<kd*cell<<','<<gd*cell<<','<<vd*cell<<','
       <<ka*cell<<','<<ga*cell<<','<<coupling*cell<<','<<volume<<','<<radius<<','
       <<axFv*cell<<','<<axWall*cell<<','<<axTv*cell<<','<<axFar*cell<<','<<(fv>0.?phi2Fv/fv:0.)<<','<<aMax;
    for(double x:centroid)out<<','<<x;
    for(double x:second)out<<','<<x;
    out<<','<<dt*frequency;
    if(monitor.enabled())for(double value:monitor.measure(v,m.lambda,m.deltaV,m.mass2,m.massPower,m.vacuumStiffening,m.formationParameters()))out<<','<<value;
    if(tracking.enabled())for(double value:tracked)out<<','<<value;
    out<<'\n';
}
int main(int argc,char** argv){try{
    SessionGuard guard(argc,argv,false);ParameterParser p(argc,argv);numeric(p);RunParameters<double> run(p);
    const int every=p.get<int>("output_every",10)(),frames=p.get<int>("snapshot_every",0)(),
              slices=p.get<int>("slice_every",0)();
    const double margin=p.get<double>("far_margin",6.)();
    if(run.N<4||run.N>384||!(margin>=0.)||run.eType!=PV2||run.t0!=0.||run.a0!=1.||run.withGWs||run.doResolutionPreserving
       ||run.expansion||run.dt<=0.||run.lSide<=0.||run.dt/run.dx>.25||run.tMax<=0.||every<=0||frames<0||slices<0)
        throw std::invalid_argument("Require N=4..384, PV2, static background, positive grid/time, dt/dx<=.25");
    PocketMonitor monitor(run.N,run.lSide,p.get<double>("monitor_radius",0.)());
    TrackingMonitor tracking(run.N,run.lSide,p.get<double>("tracking_inner",0.)(),p.get<double>("tracking_outer",0.)());
    const double stepsFloat=run.tMax/run.dt;
    if(!std::isfinite(stepsFloat)||stepsFloat>10000000.||std::abs(stepsFloat-std::round(stepsFloat))>1e-8)
        throw std::invalid_argument("Require integral bounded steps");
    const std::filesystem::path output(p.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("result_dir exists");
    PocketInitialState initial(p.get<std::string>("initial_state")(),run.N,run.lSide);
    auto toolbox=MemoryToolBox<3>::makeShared(run.N,1);PocketModel model(p,run,toolbox,initial);
    initial.release();
    std::filesystem::create_directories(output);
    PositionVerlet<double> integrator(model,run);std::ofstream history(output/"history.csv");
    history<<std::setprecision(17)<<"u,total,dilaton_kinetic,dilaton_gradient,dilaton_potential,"
           "axion_kinetic,axion_gradient,coupling,fv_volume,fv_radius,axion_fv,axion_wall,axion_tv,axion_far,"
           "phi2_fv,a_max,cx,cy,cz,sxx,syy,szz,sxy,sxz,syz,dt_frequency_bound";
    if(monitor.enabled())history<<",monitor_energy,monitor_phi_energy,monitor_axion_energy,monitor_phi_flux,monitor_axion_flux,monitor_exchange";
    if(tracking.enabled())history<<",tracking_energy,tracking_phi_energy,tracking_axion_energy,tracking_phi_flux,tracking_axion_flux,tracking_exchange,tracking_lab_phi_flux,tracking_lab_axion_flux,tracking_phi_advection,tracking_axion_advection,tracking_cx,tracking_cy,tracking_cz,tracking_vx,tracking_vy,tracking_vz,tracking_phase_volume,tracking_coherence_x,tracking_coherence_y,tracking_coherence_z,tracking_max_fv_radius";
    history<<'\n';
    const int steps=std::round(stepsFloat);
    for(int step=0;step<=steps;++step){const double u=step*run.dt;
        const bool record=step%every==0||step==steps,frame=frames>0&&step%frames==0,plane=slices>0&&step%slices==0;
        if(record||frame||plane||step==0||step==steps){
            integrator.sync(model,u);const auto v=state(model,run.N);
            if(record)measure(model,v,history,run.N,run.lSide,u,run.dt,margin,monitor,tracking);
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
