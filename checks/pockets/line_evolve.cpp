// Coupled planar fields using actual CosmoLattice NDim=1 and PV2.
// Energies are per transverse area; no spherical stability/abundance inference.
#include "TempLat/session/sessionguard.h"
#include "line_model.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include "accurate_sum.h"
#include "line_window.h"
#include "line_wall_geometry.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace TempLat;
static std::vector<double> state(PocketLineModel& model,int n) {
    auto phi=model.fldS(0_c).getFullNDHostView(),a=model.fldS(1_c).getFullNDHostView();
    auto pp=model.piS(0_c).getFullNDHostView(),pa=model.piS(1_c).getFullNDHostView();
    const auto layout=model.fldS(0_c).getLayout();std::vector<double> v(4*n);
    for(int j=0;j<n;++j) {
        device::IdxArray<1> p;layout.putMemoryIndexFromSpatialLocationInto(p,j);
        v[4*j]=phi(p[0]);v[4*j+1]=a(p[0]);v[4*j+2]=pp(p[0]);v[4*j+3]=pa(p[0]);
    }
    for(double x:v)if(!std::isfinite(x))throw std::runtime_error("Nonfinite line field");
    return v;
}
static void save(const auto& v,const std::filesystem::path& path) {
    std::ofstream file(path);file<<std::setprecision(17);
    for(std::size_t j=0;j<v.size();j+=4)file<<v[j]<<' '<<v[j+1]<<' '<<v[j+2]<<' '<<v[j+3]<<'\n';
    if(!file)throw std::runtime_error("Line state write failed");
}
int main(int argc,char** argv) {try {
    SessionGuard guard(argc,argv,false);ParameterParser parser(argc,argv);RunParameters<double> run(parser);
    const bool wallDiagnostics=parser.get<bool>("line_wall_diagnostics",false)();
    const double slabPadding=parser.get<double>("line_slab_padding",0.)();
    if(!std::isfinite(slabPadding)||slabPadding<0.)throw std::invalid_argument("Invalid slab padding");
    const int every=parser.get<int>("output_every",1)();
    const int spatialOrder=parser.get<int>("line_spatial_order",2)();
    if(spatialOrder!=2&&spatialOrder!=4)throw std::invalid_argument("line_spatial_order must be 2 or 4");
    if(const auto it=parser.getParams().find("line_spatial_order");it!=parser.getParams().end()) {
        std::istringstream text(it->second);int value;
        if(!(text>>value))throw std::invalid_argument("Invalid spatial order");
        text>>std::ws;if(!text.eof())throw std::invalid_argument("Trailing spatial order text");
    }
    LineWindow window(run.N,run.lSide,parser.get<double>("line_window_inner",0.)(),
                      parser.get<double>("line_window_outer",0.)(),
                      parser.get<double>("line_window_centre",0.)(),
                      parser.get<double>("line_window_velocity",0.)(),spatialOrder);
    const double stepsFloat=run.tMax/run.dt;
    if(run.N<4||run.N>16384||run.eType!=PV2||run.expansion||run.withGWs||run.doResolutionPreserving||
       run.t0!=0.||run.a0!=1.||!std::isfinite(run.dt)||run.dt<=0.||!std::isfinite(run.lSide)||run.lSide<=0.||
       run.dt/run.dx>.25||!std::isfinite(stepsFloat)||stepsFloat<=0.||stepsFloat>10000000.||
       std::abs(stepsFloat-std::round(stepsFloat))>1e-8||every<=0)
        throw std::invalid_argument("Require N=4..16384, static PV2 and positive integral steps, dt/dx<=.25");
    const std::filesystem::path output(parser.get<std::string>("result_dir")());
    if(std::filesystem::exists(output))throw std::invalid_argument("Line output exists");
    std::ifstream input(parser.get<std::string>("initial_state")());std::vector<double> initial(4*run.N);
    for(double& x:initial)if(!(input>>x)||!std::isfinite(x))throw std::invalid_argument("Invalid line initial data");
    input>>std::ws;if(!input.eof())throw std::invalid_argument("Trailing line initial data");
    auto toolbox=MemoryToolBox<1>::makeShared(run.N,spatialOrder==4?2:1);PocketLineModel model(parser,run,toolbox,initial);
    PositionVerlet<double> integrator(model,run);std::filesystem::create_directories(output);
    std::ofstream history(output/"history.csv");history<<std::setprecision(17)<<"t,total,phi_energy,axion_energy,exchange,dt_frequency_bound";
    if(window.enabled())history<<",window_energy,window_phi_energy,window_axion_energy,window_phi_flux,window_axion_flux,window_exchange,window_phi_lab_flux,window_axion_lab_flux,window_phi_advection,window_axion_advection";
    if(wallDiagnostics)history<<",wall_x,wall_phase_velocity,wall_width,wall_width_cells";
    if(slabPadding>0.)history<<",slab_left,slab_right,slab_left_velocity,slab_right_velocity,slab_min_width_cells,slab_energy,slab_phi_energy,slab_axion_energy,slab_phi_flux,slab_axion_flux,slab_exchange,slab_phi_lab_flux,slab_axion_lab_flux,slab_phi_advection,slab_axion_advection";
    history<<'\n';
    const int steps=std::round(stepsFloat);
    for(int step=0;step<=steps;++step) {
        if(step%every==0||step==steps) {
            integrator.sync(model,step*run.dt);const auto v=state(model,run.N);
            AccurateSum ephi,ea,exchange;double curvature=0.;
            for(int j=0;j<run.N;++j) {
                const auto p=4*j,q=4*((j+1)%run.N);
                const auto terms=pockets::formationPotential(v[p],v[p+1],model.formation);
                if(spatialOrder==2) {
                    ephi+=.5*v[p+2]*v[p+2]+terms.vacuum+.5*std::pow((v[q]-v[p])/run.dx,2);
                    ea+=.5*v[p+3]*v[p+3]+terms.axion+.5*std::pow((v[q+1]-v[p+1])/run.dx,2);
                } else {
                    ephi+=.5*v[p+2]*v[p+2]+terms.vacuum+lineFourthGradient(v,j,0,run.N,run.dx,false);
                    ea+=.5*v[p+3]*v[p+3]+terms.axion+lineFourthGradient(v,j,1,run.N,run.dx,false);
                }
                exchange+=terms.axion_exchange(v[p+2]);
                curvature=std::max(curvature,.5*(terms.hphi_phi+terms.ha_a+
                    std::hypot(terms.hphi_phi-terms.ha_a,2.*terms.hphi_a)));
            }
            const double frequency=run.dt*std::sqrt((spatialOrder==4?16./3.:4.)/(run.dx*run.dx)+curvature);
            if(frequency>.5)throw std::runtime_error("Line dt*frequency>.5");
            history<<step*run.dt<<','<<(ephi+ea)*run.dx<<','<<double(ephi)*run.dx<<','<<double(ea)*run.dx<<','<<double(exchange)*run.dx<<','<<frequency;
            if(window.enabled())for(double value:window.measure(v,step*run.dt,model.formation))history<<','<<value;
            if(wallDiagnostics)for(double value:lineWallGeometry(v,run.N,run.lSide))history<<','<<value;
            if(slabPadding>0.) {
                const auto right=lineWallGeometry(v,run.N,run.lSide);
                auto reflected=v;
                for(int j=0;j<run.N;++j)for(int c=0;c<4;++c)reflected[4*j+c]=v[4*(run.N-1-j)+c];
                const auto mirror=lineWallGeometry(reflected,run.N,run.lSide);
                const double left=run.lSide-mirror[0],lv=-mirror[1],rv=right[1];
                const double radius=(right[0]-left)/2.;
                history<<','<<left<<','<<right[0]<<','<<lv<<','<<rv<<','<<std::min(right[3],mirror[3]);
                for(double value:window.measureSlab(v,(left+right[0])/2.,radius,(lv+rv)/2.,(rv-lv)/2.,slabPadding,model.formation))history<<','<<value;
            }
            history<<'\n';
            if(step==0)save(v,output/"initial.txt");
            if(step==steps)save(v,output/"final.txt");
        }
        if(step<steps)integrator.evolve(model,step*run.dt,false);
    }
    std::ofstream parameters(output/"parameters.txt");parameters<<parser;
    if(!history||!parameters)throw std::runtime_error("Line history write failed");
    std::cout<<"Completed planar coupled fields; no 3D pocket claim\n";
} catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}}
