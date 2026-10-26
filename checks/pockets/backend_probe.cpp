#include "TempLat/session/sessionguard.h"
#include <Kokkos_Core.hpp>
#include <iostream>
#include "accurate_sum.h"

int main(int argc,char** argv) {
    TempLat::SessionGuard guard(argc,argv,false);
    Kokkos::DefaultExecutionSpace execution;
    AccurateSum sum;sum+=1e16;sum+=1.;sum+=-1e16;
    std::cout<<"POCKET_BACKEND {\"name\":\""<<execution.name()
             <<"\",\"concurrency\":"<<execution.concurrency()
             <<",\"compensated_sum\":"<<double(sum)<<"}\n";
}
