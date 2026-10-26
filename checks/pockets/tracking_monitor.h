#ifndef POCKETS_TRACKING_MONITOR_H
#define POCKETS_TRACKING_MONITOR_H
#include <array>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "accurate_sum.h"
#include "pocket_potential.h"

// Diagnostic only; TRACKING_MONITOR_CONTRACT.md derives the moving-window
// current. Compact windows avoid the periodic unwrap cut entirely.
class TrackingMonitor {
    int n_;
    double box_,dx_,inner_,outer_;
    bool enabled_;
    std::vector<double> x_,sine_,cosine_;
    std::size_t index(int i,int j,int k) const {
        return (std::size_t(i%n_)*n_+j%n_)*n_+k%n_;
    }
    double periodic(double x) const {return x-box_*std::floor(x/box_);}
    double nearest(double x) const {return x-box_*std::floor(x/box_+.5);}
public:
    TrackingMonitor(int n,double box,double inner,double outer)
      :n_(n),box_(box),dx_(box/n),inner_(inner),outer_(outer),enabled_(inner!=0.||outer!=0.) {
        if(!std::isfinite(inner)||!std::isfinite(outer)||
           (enabled_&&!(0.<inner&&inner<outer&&outer<box/2.)))
            throw std::invalid_argument("Require tracking radii both zero or 0<inner<outer<L/2");
        if(!enabled_)return;
        x_.resize(n);sine_.resize(n);cosine_.resize(n);
        for(int i=0;i<n;++i){x_[i]=(i+.5)*dx_;sine_[i]=std::sin(2.*M_PI*x_[i]/box);cosine_[i]=std::cos(2.*M_PI*x_[i]/box);}
    }
    bool enabled() const {return enabled_;}
    template<class State>
    std::array<double,21> measure(const State& v,double lambda,double dv,double mass2,int power,double stiffening,const pockets::FormationParameters* formation=nullptr) const {
        std::array<std::vector<AccurateSum>,3> weights,dots;
        for(int a=0;a<3;++a){weights[a].resize(n_);dots[a].resize(n_);}
        bool hasFalseVacuum=false;
        for(int i=0;i<n_;++i)for(int j=0;j<n_;++j)for(int k=0;k<n_;++k){
            const auto p=4*index(i,j,k);const double phi=v[p];
            hasFalseVacuum=hasFalseVacuum||phi<.5;
            const double s=std::max(0.,std::min(1.,(phi-.1)/.8));
            const double eta=1.-s*s*s*(10.+s*(-15.+6.*s));
            const double dot=-30.*s*s*(1.-s)*(1.-s)/.8*v[p+2];
            if(eta!=0.){weights[0][i]+=eta;weights[1][j]+=eta;weights[2][k]+=eta;}
            if(dot!=0.){dots[0][i]+=dot;dots[1][j]+=dot;dots[2][k]+=dot;}
        }
        if(!hasFalseVacuum)throw std::runtime_error("Tracking monitor has no false-vacuum body");
        AccurateSum totalSum;for(const auto& w:weights[0])totalSum+=double(w);
        const double total=totalSum;
        if(!(total>0.))throw std::runtime_error("Tracking monitor has zero phase volume");
        std::array<double,3> centre{},velocity{},coherence{};
        std::array<std::vector<double>,3> coordinate;
        for(int axis=0;axis<3;++axis){
            AccurateSum sine,cosine;
            for(int i=0;i<n_;++i){sine+=double(weights[axis][i])*sine_[i];cosine+=double(weights[axis][i])*cosine_[i];}
            coherence[axis]=std::hypot(double(sine),double(cosine))/total;
            if(coherence[axis]<1e-8)throw std::runtime_error("Ambiguous tracking centroid");
            const double circular=periodic(std::atan2(double(sine),double(cosine))*box_/(2.*M_PI));
            AccurateSum correction;
            for(int i=0;i<n_;++i)correction+=double(weights[axis][i])*nearest(x_[i]-circular);
            centre[axis]=periodic(circular+double(correction)/total);
            coordinate[axis].resize(n_);AccurateSum numerator;
            for(int i=0;i<n_;++i){
                const double d=nearest(x_[i]-centre[axis]);coordinate[axis][i]=d;
                if(double(weights[axis][i])>0.&&std::abs(d)>box_/2.-dx_/2.)
                    throw std::runtime_error("Tracking phase crosses its unwrap cut");
                numerator+=double(dots[axis][i])*d;
            }
            velocity[axis]=double(numerator)/total;
        }
        double maximumFalseRadius=0.;
        const auto planeSize=std::size_t(n_)*n_;
        std::vector<double> window(planeSize),windowDot(planeSize),nextWindow(planeSize),nextWindowDot(planeSize);
        auto fillPlane=[&](int i,auto& weights,auto& derivatives){
          for(int j=0;j<n_;++j)for(int k=0;k<n_;++k){
            const auto p=std::size_t(j)*n_+k;const double x=coordinate[0][i],y=coordinate[1][j],z=coordinate[2][k];
            const double r=std::sqrt(x*x+y*y+z*z),s=std::max(0.,std::min(1.,(r-inner_)/(outer_-inner_)));
            weights[p]=1.-s*s*s*(10.+s*(-15.+6.*s));
            const double radial=-30.*s*s*(1.-s)*(1.-s)/(outer_-inner_);
            derivatives[p]=r>0.?-radial*(x*velocity[0]+y*velocity[1]+z*velocity[2])/r:0.;
          }
        };
        fillPlane(0,window,windowDot);fillPlane(1%n_,nextWindow,nextWindowDot);
        AccurateSum ePhi,eA,pPhi,pA,advPhi,advA,exchange;
        for(int i=0;i<n_;++i){
          for(int j=0;j<n_;++j)for(int k=0;k<n_;++k){
            const auto cell=index(i,j,k),p=4*cell,plane=std::size_t(j)*n_+k;const double w=window[plane],wd=windowDot[plane];
            const double phi=v[p],a=v[p+1];
            if(phi<.5){
                const double x=coordinate[0][i],y=coordinate[1][j],z=coordinate[2][k];
                maximumFalseRadius=std::max(maximumFalseRadius,std::sqrt(x*x+y*y+z*z));
            }
            const std::array<std::size_t,3> plus={index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)};
            if(w!=0.||wd!=0.){
                const double phi2=phi*phi,coupling=formation?pockets::formationPotential(phi,a,*formation).axion:.5*mass2*(power==4?phi2*phi2:phi2)*a*a;
                double densityPhi=.5*v[p+2]*v[p+2]+TempLat::pocketValue<0>(phi,a,lambda,dv,mass2,power,stiffening,formation)-coupling;
                double densityA=.5*v[p+3]*v[p+3]+coupling;
                const std::array<std::size_t,3> minus={index(i+n_-1,j,k),index(i,j+n_-1,k),index(i,j,k+n_-1)};
                for(int axis=0;axis<3;++axis){
                    const auto q=4*plus[axis],r=4*minus[axis];
                    densityPhi+=((v[q]-phi)*(v[q]-phi)+(phi-v[r])*(phi-v[r]))/(4.*dx_*dx_);
                    densityA+=((v[q+1]-a)*(v[q+1]-a)+(a-v[r+1])*(a-v[r+1]))/(4.*dx_*dx_);
                }
                ePhi+=w*densityPhi;eA+=w*densityA;advPhi+=wd*densityPhi;advA+=wd*densityA;
                exchange+=formation?w*pockets::formationPotential(phi,a,*formation).axion_dphi*v[p+2]:w*mass2*(power==4?2.*phi2*phi:phi)*a*a*v[p+2];
            }
            const std::array<double,3> followingWeights={nextWindow[plane],window[std::size_t((j+1)%n_)*n_+k],window[std::size_t(j)*n_+(k+1)%n_]};
            for(int axis=0;axis<3;++axis){
                const auto next=plus[axis];const double difference=w-followingWeights[axis];
                if(difference==0.)continue;
                const auto q=4*next;const double factor=-.5*dx_*difference;
                pPhi+=factor*(v[q]-phi)*(v[q+2]+v[p+2]);
                pA+=factor*(v[q+1]-a)*(v[q+3]+v[p+3]);
            }
          }
          if(i+1<n_){
            window.swap(nextWindow);windowDot.swap(nextWindowDot);
            fillPlane((i+2)%n_,nextWindow,nextWindowDot);
          }
        }
        const double cell=dx_*dx_*dx_,phiAdvection=double(advPhi)*cell,aAdvection=double(advA)*cell;
        return {(double(ePhi)+double(eA))*cell,double(ePhi)*cell,double(eA)*cell,
            double(pPhi)-phiAdvection,double(pA)-aAdvection,double(exchange)*cell,
            double(pPhi),double(pA),phiAdvection,aAdvection,centre[0],centre[1],centre[2],
            velocity[0],velocity[1],velocity[2],total*cell,coherence[0],coherence[1],coherence[2],maximumFalseRadius};
    }
};
#endif
