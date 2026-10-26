#include "TempLat/session/sessionguard.h"
#include "bubble_model.h"
#include "../pockets/initial_state.h"
#include "../pockets/accurate_sum.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace TempLat;
#include "state_io.h"

template<class State>
static void measure(BubbleModel& m,const State& v,std::ofstream& out,int n,double box,double t,double dt) {
    const double dx=box/n,cell=dx*dx*dx;
    AccurateSum es,ea,exchange,ks,gs,us,ka,ga,va,smean,ssecond;
    double bound=0,smin=std::numeric_limits<double>::infinity(),smax=-std::numeric_limits<double>::infinity();
    std::size_t falseHalf=0,falseBarrier=0,nearTrue=0;
    const double barrier=m.lambda/(2*m.lambda+24*m.deltaV);
    auto index=[n](int i,int j,int k){return 4*((std::size_t(i%n)*n+j%n)*n+k%n);};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k) {
        const auto p=index(i,j,k); const double s=v[p],a=v[p+1],ps=v[p+2],pa=v[p+3];
        smean+=s;ssecond+=s*s;smin=std::min(smin,s);smax=std::max(smax,s);
        falseHalf+=s<.5;falseBarrier+=s<barrier;nearTrue+=std::abs(s-1)<.1;
        double gS=0,gA=0;
        for(auto z:{index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)}){
            gS+=std::pow(v[z]-s,2)/(2*dx*dx);gA+=std::pow(v[z+1]-a,2)/(2*dx*dx);
        }
        const double q=1-s,shape=2*std::pow(std::sin(a/(2*m.fa)),2);
        const double u=q*q*(m.lambda*s*s/4+m.deltaV*(1+2*s+3*s*s));
        const double w=m.fa*m.fa*m.mass*m.mass*s*s*shape;
        ks+=ps*ps/2; gs+=gS; us+=u; ka+=pa*pa/2; ga+=gA; va+=w;
        es+=ps*ps/2+gS+u;ea+=pa*pa/2+gA+w;
        exchange+=2*m.fa*m.fa*m.mass*m.mass*s*ps*shape;
        const double hss=m.lambda/2*(1-6*s+6*s*s)+m.deltaV*(36*s*s-24*s)+2*m.fa*m.fa*m.mass*m.mass*shape;
        const double haa=m.mass*m.mass*s*s*std::cos(a/m.fa);
        const double hsa=2*m.fa*m.mass*m.mass*s*std::sin(a/m.fa);
        bound=std::max(bound,.5*(hss+haa+std::hypot(hss-haa,2*hsa)));
    }
    const double frequency=dt*std::sqrt(12/(dx*dx)+bound);
    if(frequency>.5)throw std::runtime_error("dt frequency bound > .5");
    out<<t<<','<<(es+ea)*cell<<','<<es*cell<<','<<ea*cell<<','<<exchange*cell
       <<','<<ks*cell<<','<<gs*cell<<','<<us*cell<<','<<ka*cell<<','<<ga*cell<<','<<va*cell<<','<<frequency;
    const double count=double(n)*n*n,mean=double(smean)/count;
    out<<','<<smin<<','<<smax<<','<<mean<<','<<std::sqrt(std::max(0.,double(ssecond)/count-mean*mean))
       <<','<<falseHalf/count<<','<<falseBarrier/count<<','<<nearTrue/count<<'\n';
}
int main(int argc,char** argv){try {
    SessionGuard guard(argc,argv,false);ParameterParser p(argc,argv);RunParameters<double> run(p);
    // Reject malformed numerics before parser conversions and file allocation.
    for(const char* key:{"N","lSide","dt","tMax","lambda","delta_v","mass_tv","axion_fa","output_every","snapshot_every"}) {
        const auto it=p.getParams().find(key);if(it==p.getParams().end())throw std::invalid_argument("Missing numeric parameter");
        std::istringstream in(it->second);double x;
        if(!(in>>x)||!std::isfinite(x))throw std::invalid_argument("Invalid numeric parameter");
        in>>std::ws;if(!in.eof())throw std::invalid_argument("Trailing numeric parameter text");
        if((std::string(key)=="N"||std::string(key)=="output_every"||std::string(key)=="snapshot_every")&&(x!=std::floor(x)||x<0||x>10000000))throw std::invalid_argument("Invalid integer parameter");
    }
    const int every=p.get<int>("output_every")(),frames=p.get<int>("snapshot_every")();
    const double stepsFloat=run.tMax/run.dt;
    if(run.N<4||run.N>192||run.eType!=PV2||run.t0!=0||run.a0!=1||run.withGWs||run.expansion||run.doResolutionPreserving||
       run.dt<=0||run.lSide<=0||run.dt/run.dx>.25||run.tMax<=0||every<1||frames<1||!std::isfinite(stepsFloat)||
       stepsFloat>1000000||std::abs(stepsFloat-std::round(stepsFloat))>1e-8)throw std::invalid_argument("Invalid single-bubble run contract");
    const std::filesystem::path output(p.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("Output exists");
    PocketInitialState initial(p.get<std::string>("initial_state")(),run.N,run.lSide);
    auto toolbox=MemoryToolBox<3>::makeShared(run.N,1);BubbleModel model(p,run,toolbox,initial);initial.release();
    std::filesystem::create_directories(output);
    PositionVerlet<double> integrator(model,run);
    std::ofstream h(output/"history.csv");h<<std::setprecision(17)<<"t,total,E_S,E_A,Q_A,K_S,G_S,U,K_A,G_A,V_A,dt_frequency,S_min,S_max,S_mean,S_std,F_half,F_barrier,near_true_fraction\n";
    const int steps=std::round(stepsFloat);
    for(int step=0;step<=steps;++step){const double t=step*run.dt;
        if(step%every==0||step%frames==0||step==steps){
            integrator.sync(model,t);const auto v=state(model,run.N);
            if(step%every==0||step==steps)measure(model,v,h,run.N,run.lSide,t,run.dt);
            if(step==0)snapshot(v,output/"initial.bin",run.N,run.lSide,t);
            if(step==steps)snapshot(v,output/"final.bin",run.N,run.lSide,t);
            if(step%frames==0||step==steps){std::ostringstream name;name<<"frame_"<<std::setw(7)<<std::setfill('0')<<step<<".bin";
                if(run.N<=4)snapshot(v,output/name.str(),run.N,run.lSide,t);
                else slice(v,output/name.str(),run.N,run.lSide,t);}
        }
        if(step<steps)integrator.evolve(model,t,false);
    }
    std::ofstream pars(output/"parameters.txt");pars<<p;if(!h||!pars)throw std::runtime_error("Output failed");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
