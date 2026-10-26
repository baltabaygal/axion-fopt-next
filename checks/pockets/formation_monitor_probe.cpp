// Actual monitor headers, without a framework build or lattice evolution.
#include "pocket_monitor.h"
#include "tracking_monitor.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
int main(int argc,char** argv) {
    if(argc!=3)return 2;
    const int n=std::atoi(argv[1]),mode=std::atoi(argv[2]);
    if(n<4 || n>24 || mode<0 || mode>2)return 2;
    std::vector<double> state(4*std::size_t(n)*n*n);
    for(double& value:state)if(!(std::cin>>value))return 3;
    pockets::FormationParameters formation;
    const auto* params=mode==0?&formation:nullptr;
    const int power=mode==2?4:2;
    const double stiffening=mode==2?12.:0.;
    PocketMonitor fixed(n,8.,2.7);
    TrackingMonitor moving(n,8.,2.4,3.5);
    std::cout<<std::setprecision(17);
    for(double value:fixed.measure(state,8.,.01,1.,power,stiffening,params))std::cout<<value<<' ';
    std::cout<<'\n';
    for(double value:moving.measure(state,8.,.01,1.,power,stiffening,params))std::cout<<value<<' ';
    std::cout<<'\n';
}
