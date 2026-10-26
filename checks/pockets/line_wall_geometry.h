#ifndef POCKETS_LINE_WALL_GEOMETRY_H
#define POCKETS_LINE_WALL_GEOMETRY_H
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
// Unique rising phase interface; ambiguity is an error, never silently filtered.
inline std::array<double,4> lineWallGeometry(const std::vector<double>& state,int n,double box) {
    if(n<4||state.size()!=4*std::size_t(n)||!std::isfinite(box)||box<=0.)
        throw std::invalid_argument("Invalid line wall geometry");
    const double dx=box/n;
    const auto periodic=[box](double x){return x-box*std::floor(x/box);};
    const auto nearest=[box](double x){return x-box*std::floor(x/box+.5);};
    const auto crossing=[&](double level) {
        int from=-1,count=0;
        for(int j=0;j<n;++j)if(state[4*j]<=level&&state[4*((j+1)%n)]>level){from=j;++count;}
        if(count!=1)throw std::runtime_error("Ambiguous rising line wall crossing");
        const int to=(from+1)%n;
        const double delta=state[4*to]-state[4*from],fraction=(level-state[4*from])/delta;
        const double momentum=(1.-fraction)*state[4*from+2]+fraction*state[4*to+2];
        return std::array<double,2>{periodic((from+.5+fraction)*dx),-momentum/(delta/dx)};
    };
    const auto centre=crossing(.5),lo=crossing(.1),hi=crossing(.9);
    const double lower=nearest(lo[0]-centre[0]),upper=nearest(hi[0]-centre[0]);
    if(!(lower<0.&&upper>0.))throw std::runtime_error("Invalid line wall width ordering");
    const double width=upper-lower;
    return {centre[0],centre[1],width,width/dx};
}
#endif
