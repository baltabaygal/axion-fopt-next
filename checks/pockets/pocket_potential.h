#ifndef POCKETS_POTENTIAL_H
#define POCKETS_POTENTIAL_H
#include <cmath>
#include "formation_potential.h"
#ifdef DEVICE_INLINE_FUNCTION
#define POCKETS_SCALAR_INLINE DEVICE_INLINE_FUNCTION
#else
#define POCKETS_SCALAR_INLINE inline
#endif
namespace TempLat {
// Dilaton vev v=1: false vacuum phi=0, true vacuum phi=1 with V(1)=0, V(0)=deltaV.
// Axion mass m_a^2=mass2*phi^2. Order 0: V; 1: dV/dphi; 2: dV/da.
template<int Order>
POCKETS_SCALAR_INLINE double pocketValue(double phi, double a, double lambda, double deltaV, double mass2, int power=2,double stiffening=0.,const pockets::FormationParameters* formation=nullptr) {
    if(formation) {
        const auto v=pockets::formationPotential(phi,a,*formation);
        if constexpr(Order==0)return v.potential();
        if constexpr(Order==1)return v.dphi;
        if constexpr(Order==2)return v.da;
    }
    const double q=1.-phi;
    if(stiffening!=0.) {
        const double phi2=phi*phi,y=2.*phi-1.,y2=y*y,y4=y2*y2,y8=y4*y4;
        const double massShape=power==4?phi2*phi2:phi2;
        if constexpr(Order==0) return q*q*(.25*lambda*phi2*(1.+stiffening*y8)+deltaV*(1.+2.*phi+3.*phi2))+.5*mass2*massShape*a*a;
        if constexpr(Order==1) return .5*lambda*phi*q*(1.-2.*phi)*(1.+stiffening*y8)
            +4.*lambda*stiffening*phi2*q*q*y4*y2*y+12.*deltaV*phi2*(phi-1.)
            +mass2*a*a*(power==4?2.*phi2*phi:phi);
        if constexpr(Order==2) return mass2*massShape*a;
    }
    if(power==4) {
        const double phi2=phi*phi;
        if constexpr (Order==0) return q*q*(.25*lambda*phi2+deltaV*(1.+2.*phi+3.*phi2))+.5*mass2*phi2*phi2*a*a;
        if constexpr (Order==1) return .5*lambda*phi*q*(1.-2.*phi)+12.*deltaV*phi2*(phi-1.)+2.*mass2*phi2*phi*a*a;
        if constexpr (Order==2) return mass2*phi2*phi2*a;
    }
    if constexpr (Order==0) return q*q*(.25*lambda*phi*phi+deltaV*(1.+2.*phi+3.*phi*phi))+.5*mass2*phi*phi*a*a;
    if constexpr (Order==1) return .5*lambda*phi*q*(1.-2.*phi)+12.*deltaV*phi*phi*(phi-1.)+mass2*phi*a*a;
    if constexpr (Order==2) return mass2*phi*phi*a;
}
POCKETS_SCALAR_INLINE double pocketDilatonCurvature(double phi,double lambda,double deltaV,double stiffening) {
    const double phi2=phi*phi;
    const double base=.5*lambda*(1.-6.*phi+6.*phi2),tilt=deltaV*(36.*phi2-24.*phi);
    if(stiffening==0.)return base+tilt;
    const double q=1.-phi,y=2.*phi-1.,y2=y*y,y4=y2*y2,y6=y4*y2,y7=y6*y,y8=y4*y4;
    return base*(1.+stiffening*y8)
        +2.*(.5*lambda*phi*q*(1.-2.*phi))*(16.*stiffening*y7)
        +(.25*lambda*phi2*q*q)*(224.*stiffening*y6)+tilt;
}
}
#undef POCKETS_SCALAR_INLINE
#endif
