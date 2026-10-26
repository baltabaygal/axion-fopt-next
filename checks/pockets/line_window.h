#ifndef POCKETS_LINE_WINDOW_H
#define POCKETS_LINE_WINDOW_H
#include <algorithm>
#include <array>
#include <vector>
#include "formation_potential.h"
#include "accurate_sum.h"
#include "line_spatial_energy.h"
// Passive compact window with specified motion; never a force or wall trajectory.
class LineWindow {
    int n_,order_;double box_,dx_,inner_,outer_,centre_,velocity_;
public:
    LineWindow(int n,double box,double inner,double outer,double centre,double velocity,int order=2)
      :n_(n),order_(order),box_(box),dx_(box/n),inner_(inner),outer_(outer),centre_(centre),velocity_(velocity) {
        if((order!=2&&order!=4)||n<4||!std::isfinite(box)||box<=0.||!std::isfinite(inner)||!std::isfinite(outer)||
           !std::isfinite(centre)||!std::isfinite(velocity)||inner<0.||outer<0.||
           (outer!=0.&&!(inner<outer&&outer<box/2.))||(outer==0.&&inner!=0.))
            throw std::invalid_argument("Invalid passive line window");
    }
    bool enabled() const {return outer_>0.;}
    std::array<double,10> measure(const std::vector<double>& state,double t,
                                  const pockets::FormationParameters& p) const {
        if(!enabled()||state.size()!=4*std::size_t(n_)||!std::isfinite(t))
            throw std::invalid_argument("Invalid line window measurement");
        std::vector<double> weights(n_),dots(n_);
        for(int j=0;j<n_;++j) {
            double d=(j+.5)*dx_-centre_-velocity_*t;d-=box_*std::floor(d/box_+.5);
            const double s=std::clamp((std::abs(d)-inner_)/(outer_-inner_),0.,1.);
            weights[j]=1.-s*s*s*(10.+s*(-15.+6.*s));
            const double radial=-30.*s*s*(1.-s)*(1.-s)/(outer_-inner_);
            dots[j]=-velocity_*(d>0.?1.:d<0.?-1.:0.)*radial;
        }
        return measureWeights(state,weights,dots,p);
    }
    std::array<double,10> measureSlab(const std::vector<double>& state,
            double centre,double radius,double centreVelocity,double radiusVelocity,
            double padding,const pockets::FormationParameters& p) const {
        if(!std::isfinite(centre)||!std::isfinite(radius)||!std::isfinite(centreVelocity)||
           !std::isfinite(radiusVelocity)||!std::isfinite(padding)||
           !(padding>0.&&padding<radius&&radius+padding<box_/2.))
            throw std::invalid_argument("Invalid wall-following slab window");
        std::vector<double> weights(n_),dots(n_);
        for(int j=0;j<n_;++j) {
            double d=(j+.5)*dx_-centre;d-=box_*std::floor(d/box_+.5);
            const double z=std::clamp((std::abs(d)-radius+padding)/(2.*padding),0.,1.);
            weights[j]=1.-z*z*z*(10.+z*(-15.+6.*z));
            dots[j]=30.*z*z*(1.-z)*(1.-z)*
                ((d>0.?1.:d<0.?-1.:0.)*centreVelocity+radiusVelocity)/(2.*padding);
        }
        return measureWeights(state,weights,dots,p);
    }
private:
    std::array<double,10> measureWeights(const std::vector<double>& state,
            const std::vector<double>& weights,const std::vector<double>& dots,
            const pockets::FormationParameters& p) const {
        if(state.size()!=4*std::size_t(n_))throw std::invalid_argument("Invalid window state");
        AccurateSum ep,ea,fp,fa,ap,aa,q;
        for(int j=0;j<n_;++j) {
            const auto at=4*j,next=4*((j+1)%n_),prev=4*((j+n_-1)%n_);
            const auto v=pockets::formationPotential(state[at],state[at+1],p);
            const double gp=state[next]-state[at],gm=state[at]-state[prev];
            const double bp=state[next+1]-state[at+1],bm=state[at+1]-state[prev+1];
            double rhoP=.5*state[at+2]*state[at+2]+v.vacuum+(gp*gp+gm*gm)/(4.*dx_*dx_);
            double rhoA=.5*state[at+3]*state[at+3]+v.axion+(bp*bp+bm*bm)/(4.*dx_*dx_);
            if(order_==4) {
                rhoP=.5*state[at+2]*state[at+2]+v.vacuum+lineFourthGradient(state,j,0,n_,dx_,true);
                rhoA=.5*state[at+3]*state[at+3]+v.axion+lineFourthGradient(state,j,1,n_,dx_,true);
            }
            ep+=weights[j]*rhoP;ea+=weights[j]*rhoA;ap+=dots[j]*rhoP;aa+=dots[j]*rhoA;
            const double dw=weights[j]-weights[(j+1)%n_];
            if(order_==2) {
                fp+=-gp*(state[at+2]+state[next+2])/(2.*dx_)*dw;
                fa+=-bp*(state[at+3]+state[next+3])/(2.*dx_)*dw;
            } else for(const auto& [hop,c]:lineFourthBonds) {
                const int right=(j+hop)%n_,to=4*right;
                const double difference=weights[j]-weights[right];
                fp+=-c*(state[to]-state[at])*(state[at+2]+state[to+2])/(2.*dx_)*difference;
                fa+=-c*(state[to+1]-state[at+1])*(state[at+3]+state[to+3])/(2.*dx_)*difference;
            }
            q+=weights[j]*v.axion_exchange(state[at+2]);
        }
        const double advP=double(ap)*dx_,advA=double(aa)*dx_;
        return {(ep+ea)*dx_,double(ep)*dx_,double(ea)*dx_,double(fp)-advP,double(fa)-advA,
                double(q)*dx_,double(fp),double(fa),advP,advA};
    }
};
#endif
