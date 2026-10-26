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
#include "bubble_insertion.h"

template<class State>
static std::array<double,3> measure(BubbleModel& m,const State& v,std::ostream& out,int n,double box,double t,double dt,double h0,const std::array<double,3>& inserted) {
    const double dx=box/n,cell=dx*dx*dx,scale=m.aI;
    const double a2=scale*scale,a6=a2*a2*a2;
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
            gS+=std::pow(v[z]-s,2)/(2*dx*dx*a2);gA+=std::pow(v[z+1]-a,2)/(2*dx*dx*a2);
        }
        const double q=1-s,shape=2*std::pow(std::sin(a/(2*m.fa)),2);
        const double u=q*q*(m.lambda*s*s/4+m.deltaV*(1+2*s+3*s*s));
        const double w=m.fa*m.fa*m.mass*m.mass*s*s*shape;
        ks+=ps*ps/(2*a6); gs+=gS; us+=u; ka+=pa*pa/(2*a6); ga+=gA; va+=w;
        es+=ps*ps/(2*a6)+gS+u;ea+=pa*pa/(2*a6)+gA+w;
        exchange+=2*m.fa*m.fa*m.mass*m.mass*s*ps*shape/a2;
        const double hss=m.lambda/2*(1-6*s+6*s*s)+m.deltaV*(36*s*s-24*s)+2*m.fa*m.fa*m.mass*m.mass*shape;
        const double haa=m.mass*m.mass*s*s*std::cos(a/m.fa);
        const double hsa=2*m.fa*m.mass*m.mass*s*std::sin(a/m.fa);
        bound=std::max(bound,.5*(hss+haa+std::hypot(hss-haa,2*hsa)));
    }
    const double frequency=dt*std::sqrt(12/(dx*dx)+a2*bound);
    if(frequency>.5)throw std::runtime_error("dt frequency bound > .5");
    out<<t<<','<<(es+ea)*cell<<','<<es*cell<<','<<ea*cell<<','<<exchange*cell
       <<','<<ks*cell<<','<<gs*cell<<','<<us*cell<<','<<ka*cell<<','<<ga*cell<<','<<va*cell<<','<<frequency;
    const double count=double(n)*n*n,mean=double(smean)/count;
    out<<','<<smin<<','<<smax<<','<<mean<<','<<std::sqrt(std::max(0.,double(ssecond)/count-mean*mean))
       <<','<<falseHalf/count<<','<<falseBarrier/count<<','<<nearTrue/count
       <<','<<scale<<','<<t+h0*t*t/2<<','<<h0/scale<<','<<h0/(scale*scale)<<','<<inserted[0]<<','<<inserted[1]<<','<<inserted[2]<<'\n';
    return {(es+ea)*cell,es*cell,ea*cell};
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
    if(run.expansion && (!run.fixedBackground||std::abs(run.omegaEoS-1./3.)>1e-12||!std::isfinite(run.H0)||run.H0<=0))
        throw std::invalid_argument("Require prescribed radiation expansion with positive H0");
    const double h0=run.expansion?run.H0:0.;
    const int every=p.get<int>("output_every")(),frames=p.get<int>("snapshot_every")();
    const double stepsFloat=run.tMax/run.dt;
    if(run.N<4||run.N>192||run.eType!=PV2||run.t0!=0||run.a0!=1||run.withGWs||run.doResolutionPreserving||
       run.dt<=0||run.lSide<=0||run.dt/run.dx>.25||run.tMax<=0||every<1||frames<1||!std::isfinite(stepsFloat)||
       stepsFloat>1000000||std::abs(stepsFloat-std::round(stepsFloat))>1e-8)throw std::invalid_argument("Invalid single-bubble run contract");
    const std::filesystem::path output(p.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("Output exists");
    PocketInitialState initial(p.get<std::string>("initial_state")(),run.N,run.lSide);
    auto toolbox=MemoryToolBox<3>::makeShared(run.N,1);BubbleModel model(p,run,toolbox,initial);initial.release();
    std::filesystem::create_directories(output);
    PositionVerlet<double> integrator(model,run);
    std::ofstream h(output/"history.csv");h<<std::setprecision(17)<<"u,total,E_S,E_A,Q_A,K_S,G_S,U,K_A,G_A,V_A,dt_frequency,S_min,S_max,S_mean,S_std,F_half,F_barrier,near_true_fraction,a,t_cosmic,H_conformal,H_physical,I_total,I_S,I_A\n";
    const auto events=readBirths(p.get<std::string>("candidate_file", "none")(),run.lSide,run.tMax);
    const double radius=p.get<double>("seed_radius",8.)(),width=p.get<double>("seed_width",2.)(),tolerance=p.get<double>("core_tolerance",.05)();
    if(!std::isfinite(radius)||!std::isfinite(width)||!std::isfinite(tolerance)||radius<=0||width<=0||tolerance<0)
        throw std::invalid_argument("Invalid insertion profile");
    if(!events.empty()&&radius+4*width>run.lSide/2)throw std::invalid_argument("Initial seed support exceeds periodic half box");
    std::ofstream ledger(output/"births.csv");ledger<<std::setprecision(17)<<"u,x,y,z,accepted,max_core_abs_S,core_cells,delta_total,delta_S,delta_A\n";
    std::array<double,3> inserted{0,0,0};std::size_t eventIndex=0;int historyIndex=0,frameIndex=0,steps=0;
    double u=0;const double baseDt=run.dt,cadence=every*baseDt,frameCadence=frames*baseDt,epsilon=1e-11;
    for(;;){
        integrator.sync(model,u);
        while(eventIndex<events.size()&&events[eventIndex].u<=u+epsilon){
            const auto& event=events[eventIndex++];const double a=1+h0*u;const auto current=state(model,run.N);
            const auto check=inspectCore(current,run.N,run.lSide,a,event,radius,tolerance);
            std::array<double,3> delta{0,0,0};
            if(check.accepted){
                const auto before=measure(model,current,h,run.N,run.lSide,u,baseDt,h0,inserted);
                snapshot(current,output/("birth_pre_"+std::to_string(eventIndex)+".bin"),run.N,run.lSide,u);
                insertBirth(model,run.N,run.lSide,a,event,radius,width);
                const auto afterState=state(model,run.N);std::ostringstream sink;
                const auto after=measure(model,afterState,sink,run.N,run.lSide,u,baseDt,h0,inserted);
                for(int d=0;d<3;++d){delta[d]=after[d]-before[d];inserted[d]+=delta[d];}
                measure(model,afterState,h,run.N,run.lSide,u,baseDt,h0,inserted);
                snapshot(afterState,output/("birth_post_"+std::to_string(eventIndex)+".bin"),run.N,run.lSide,u);
            }
            ledger<<u<<','<<event.x<<','<<event.y<<','<<event.z<<','<<check.accepted<<','<<check.maxCore<<','<<check.cells;
            for(double d:delta)ledger<<','<<d;ledger<<'\n';
        }
        const bool record=u+epsilon>=historyIndex*cadence,frame=u+epsilon>=frameIndex*frameCadence,final=u+epsilon>=run.tMax;
        if(record||frame||final||u==0){const auto v=state(model,run.N);
            if(record||final){measure(model,v,h,run.N,run.lSide,u,baseDt,h0,inserted);if(record)++historyIndex;}
            if(u==0)snapshot(v,output/"initial.bin",run.N,run.lSide,u);
            if(final)snapshot(v,output/"final.bin",run.N,run.lSide,u);
            if(frame||final){std::ostringstream name;name<<"frame_"<<std::setw(7)<<std::setfill('0')<<frameIndex<<".bin";
                if(run.N<=4)snapshot(v,output/name.str(),run.N,run.lSide,u);else slice(v,output/name.str(),run.N,run.lSide,u);
                if(frame)++frameIndex;}
        }
        if(final)break;
        double target=std::min({u+baseDt,run.tMax,historyIndex*cadence,frameIndex*frameCadence});
        if(eventIndex<events.size())target=std::min(target,events[eventIndex].u);
        if(target-u<=epsilon||++steps>2000000)throw std::runtime_error("Invalid event-aligned step");
        model.dt=target-u;integrator.evolve(model,u,false);u=target;
    }
    std::ofstream pars(output/"parameters.txt");pars<<p;if(!h||!pars||!ledger)throw std::runtime_error("Output failed");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
