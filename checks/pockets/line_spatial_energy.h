#ifndef POCKETS_LINE_SPATIAL_ENERGY_H
#define POCKETS_LINE_SPATIAL_ENERGY_H
#include <array>
#include <utility>
#include <vector>
inline constexpr std::array<std::pair<int,double>,2> lineFourthBonds={{{1,4./3.},{2,-1./12.}}};
inline double lineFourthGradient(const std::vector<double>& v,int j,int component,int n,double dx,bool shared) {
    double value=0.;const double f=v[4*j+component];
    for(const auto& [hop,c]:lineFourthBonds) {
        const double plus=v[4*((j+hop)%n)+component]-f;
        if(shared) {
            const double minus=f-v[4*((j+n-hop)%n)+component];
            value+=c*(plus*plus+minus*minus)/(4.*dx*dx);
        } else value+=c*plus*plus/(2.*dx*dx);
    }
    return value;
}
#endif
