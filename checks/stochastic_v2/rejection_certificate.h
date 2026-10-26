#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

// Conservative displacement envelope, not an extrapolated field prediction.
// Bootstrap |S|<=max|S0|+1. A candidate is skipped only if the resulting
// |delta S| bound is <1 and an actual sampled cell in its core stays |S|>tol.
template<class State,class Model>
double displacementEnvelope(const State& v,const Model& m,double dx,double a0,double amax,double du) {
    double s0=0,p0=0;
    for(std::size_t i=0;i<v.size();i+=4){s0=std::max(s0,std::abs(v[i]));p0=std::max(p0,std::abs(v[i+2]));}
    const double b=s0+1;
    const double force=m.lambda/4*(2*b+6*b*b+4*b*b*b)+m.deltaV*(12*b*b+12*b*b*b)+4*m.fa*m.fa*m.mass*m.mass*b;
    const double dp=amax*amax*12*b/(dx*dx)+std::pow(amax,4)*force;
    // Includes the full kick, hence also bounds the PV2 half kick. Substeps
    // without insertion obey the same sum-of-forces envelope by induction.
    return du*(p0+du*dp)/(a0*a0);
}
template<class State>
double rejectionWitness(const State& v,int n,double box,double amax,const BirthEvent& e,double radius,double tolerance,double envelope) {
    if(!std::isfinite(envelope)||envelope>=1)return -1;
    const double dx=box/n,offset=radius/(2*amax);const double c[3]={e.x,e.y,e.z};
    for(int ox:{0,-1,1})for(int oy:{0,-1,1})for(int oz:{0,-1,1}) {
        const int offsets[3]={ox,oy,oz};int ix[3];double d2=0;
        for(int d=0;d<3;++d){const int raw=int(std::floor((c[d]+offsets[d]*offset)/dx));ix[d]=(raw%n+n)%n;
            double distance=std::abs((ix[d]+.5)*dx-c[d]);distance=std::min(distance,box-distance);d2+=distance*distance;}
        if(amax*amax*d2>radius*radius)continue;
        const double lower=std::abs(v[4*((std::size_t(ix[0])*n+ix[1])*n+ix[2])])-envelope;
        // Roundoff margin deliberately favors the ordinary event path.
        if(lower>tolerance+1e-10)return lower;
    }
    return -1;
}

struct EmptyCertificate {bool valid;int block;double minimum;};
template<class State>
EmptyCertificate certifyEmpty(const State& v,int n,double box,double amax,double radius,double tolerance,double envelope) {
    if(!std::isfinite(envelope)||envelope>=1)return {false,0,0};
    const double dx=box/n;
    // Each centre belongs to the tile of its nearest cell. Its distance to
    // ANY witness cell in that tile is <=sqrt(3)*(block-.5)*dx.
    const int block=std::min(n,int(std::floor(radius/(amax*std::sqrt(3.)*dx)+.5-1e-12)));
    if(block<1)return {false,0,0};
    const int side=(n+block-1)/block;std::vector<double> maximum(std::size_t(side)*side*side,0.);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
        const auto tile=(std::size_t(i/block)*side+j/block)*side+k/block;
        maximum[tile]=std::max(maximum[tile],std::abs(v[4*((std::size_t(i)*n+j)*n+k)]));
    }
    double minimum=std::numeric_limits<double>::infinity();
    for(double value:maximum)minimum=std::min(minimum,value-envelope);
    return {minimum>tolerance+1e-10,block,minimum};
}
