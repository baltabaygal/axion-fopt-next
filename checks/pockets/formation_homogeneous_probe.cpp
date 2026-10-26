// Standalone homogeneous scalar-kernel check, not CosmoLattice evidence.
#include "formation_potential.h"
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
using State=std::array<double,5>; // phi, a, phi_dot, a_dot, integrated Q_a
State rhs(const State& y) {
    const auto v=pockets::formationPotential(y[0],y[1],{});
    return {y[2],y[3],-v.dphi,-v.da,v.axion_exchange(y[2])};
}
State shifted(const State& y,const State& k,double h) {
    State out{};for(int j=0;j<5;++j)out[j]=y[j]+h*k[j];return out;
}
void emit(double t,const State& y) {
    const auto v=pockets::formationPotential(y[0],y[1],{});
    const double ea=.5*y[3]*y[3]+v.axion;
    std::cout<<t;for(double x:y)std::cout<<','<<x;
    std::cout<<','<<ea+.5*y[2]*y[2]+v.vacuum<<','<<ea<<'\n';
}
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    const int steps=std::atoi(argv[1]);
    if(steps<20 || steps%20)return 2;
    const double h=4./steps;
    // Finite false-vacuum backreaction with initially stationary dilaton.
    State y{0.,2.,0.,.15,0.};
    std::cout<<std::setprecision(17)<<"t,phi,a,phi_dot,a_dot,integrated_Q,energy,axion_energy\n";
    emit(0.,y);
    for(int i=1;i<=steps;++i) {
        const auto k1=rhs(y),k2=rhs(shifted(y,k1,h/2.));
        const auto k3=rhs(shifted(y,k2,h/2.)),k4=rhs(shifted(y,k3,h));
        for(int j=0;j<5;++j)y[j]+=h/6.*(k1[j]+2.*k2[j]+2.*k3[j]+k4[j]);
        if(i%(steps/20)==0)emit(i*h,y);
    }
}
