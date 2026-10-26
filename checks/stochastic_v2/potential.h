#pragma once
#include <cmath>
namespace TempLat {
// Canonical S and A; theta=A/f. Full reciprocal coupling, no external mass.
template<int Order>
DEVICE_INLINE_FUNCTION double bubbleValue(double s,double a,double l,double d,double m,double f) {
    const double q=1-s, sn=std::sin(a/(2*f));
    const double shape=2*sn*sn, chi=f*f*m*m;
    if constexpr(Order==0)return q*q*(l*s*s/4+d*(1+2*s+3*s*s))+chi*s*s*shape;
    if constexpr(Order==1)return l*s*q*(1-2*s)/2-12*d*s*s*q+2*chi*s*shape;
    if constexpr(Order==2)return f*m*m*s*s*std::sin(a/f);
}
}
