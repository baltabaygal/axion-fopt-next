#ifndef POCKETS_MONITOR_H
#define POCKETS_MONITOR_H
#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "accurate_sum.h"
#include "pocket_potential.h"

// Static diagnostic only. Symmetric bond energy and its exact semi-discrete
// current; see MONITOR_CONTRACT.md. No force or boundary modification.
class PocketMonitor {
    struct Face { std::size_t from,to; double sign; };
    int n_;
    double dx_;
    bool enabled_;
    std::vector<std::size_t> inside_;
    std::vector<Face> faces_;
    std::size_t index(int i,int j,int k) const {
        return (std::size_t(i%n_)*n_+j%n_)*n_+k%n_;
    }
public:
    PocketMonitor(int n,double box,double radius):n_(n),dx_(box/n),enabled_(radius>0.) {
        if(!std::isfinite(radius)||radius<0.||radius>=box/2.)
            throw std::invalid_argument("Require 0<=monitor_radius<L/2");
        if(!enabled_)return;
        std::vector<unsigned char> mask(std::size_t(n)*n*n);
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
            const double x=(i+.5)*dx_-box/2.,y=(j+.5)*dx_-box/2.,z=(k+.5)*dx_-box/2.;
            const auto p=index(i,j,k);
            mask[p]=x*x+y*y+z*z<radius*radius;
            if(mask[p])inside_.push_back(p);
        }
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
            const auto p=index(i,j,k);
            for(auto q:{index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)})
                if(mask[p]!=mask[q])faces_.push_back({4*p,4*q,mask[p]?1.:-1.});
        }
    }
    bool enabled() const {return enabled_;}
    template<class State>
    std::array<double,6> measure(const State& v,double lambda,double dv,double mass2,int power,double stiffening=0.,const pockets::FormationParameters* formation=nullptr) const {
        AccurateSum ePhi,eA,pPhi,pA,exchange;
        const double cell=dx_*dx_*dx_;
        for(auto cellIndex:inside_){
            const auto p=4*cellIndex;
            const int i=cellIndex/(n_*n_),j=(cellIndex/n_)%n_,k=cellIndex%n_;
            const double phi=v[p],a=v[p+1],phi2=phi*phi;
            const double coupling=formation?pockets::formationPotential(phi,a,*formation).axion:.5*mass2*(power==4?phi2*phi2:phi2)*a*a;
            ePhi+=.5*v[p+2]*v[p+2]+TempLat::pocketValue<0>(phi,a,lambda,dv,mass2,power,stiffening,formation)-coupling;
            eA+=.5*v[p+3]*v[p+3]+coupling;
            const std::array<std::size_t,3> plus={index(i+1,j,k),index(i,j+1,k),index(i,j,k+1)};
            const std::array<std::size_t,3> minus={index(i+n_-1,j,k),index(i,j+n_-1,k),index(i,j,k+n_-1)};
            for(int axis=0;axis<3;++axis){
                const auto q=4*plus[axis],r=4*minus[axis];
                ePhi+=((v[q]-phi)*(v[q]-phi)+(phi-v[r])*(phi-v[r]))/(4.*dx_*dx_);
                eA+=((v[q+1]-a)*(v[q+1]-a)+(a-v[r+1])*(a-v[r+1]))/(4.*dx_*dx_);
            }
            exchange+=(formation?pockets::formationPotential(phi,a,*formation).axion_dphi:mass2*(power==4?2.*phi2*phi:phi)*a*a)*v[p+2];
        }
        for(const auto& face:faces_){
            const auto p=face.from,q=face.to;
            const double factor=-.5*dx_*face.sign;
            pPhi+=factor*(v[q]-v[p])*(v[q+2]+v[p+2]);
            pA+=factor*(v[q+1]-v[p+1])*(v[q+3]+v[p+3]);
        }
        return {(ePhi+eA)*cell,ePhi*cell,eA*cell,pPhi,pA,exchange*cell};
    }
};
#endif
