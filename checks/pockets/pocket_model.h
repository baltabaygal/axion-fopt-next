#ifndef POCKETS_MODEL_H
#define POCKETS_MODEL_H
#include "TempLat/util/staticif.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "CosmoInterface/abstractmodel.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include <cmath>
#include <vector>

namespace TempLat {
// Dilaton vev v=1: false vacuum phi=0, true vacuum phi=1 with V(1)=0, V(0)=deltaV.
// Axion mass m_a^2=mass2*phi^2. Order 0: V; 1: dV/dphi; 2: dV/da.
template<int Order>
DEVICE_INLINE_FUNCTION double pocketValue(double phi, double a, double lambda, double deltaV, double mass2) {
    const double q=1.-phi;
    if constexpr (Order==0) return q*q*(.25*lambda*phi*phi+deltaV*(1.+2.*phi+3.*phi*phi))+.5*mass2*phi*phi*a*a;
    if constexpr (Order==1) return .5*lambda*phi*q*(1.-2.*phi)+12.*deltaV*phi*phi*(phi-1.)+mass2*phi*a*a;
    if constexpr (Order==2) return mass2*phi*phi*a;
}
template<class X, class Y, int Order>
class PocketFunction : public BinaryOperator<X,Y> {
    double lambda_,deltaV_,mass2_;
public:
    PocketFunction(const X& x,const Y& y,double lambda,double deltaV,double mass2)
      : BinaryOperator<X,Y>(x,y),lambda_(lambda),deltaV_(deltaV),mass2_(mass2) {}
    template<typename... IDX>
      requires requires(std::decay_t<X> x,std::decay_t<Y> y,IDX... idx) {
        requires IsVariadicIndex<IDX...>; DoEval::eval(x,idx...); DoEval::eval(y,idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        return pocketValue<Order>(DoEval::eval(this->mR,idx...),DoEval::eval(this->mT,idx...),lambda_,deltaV_,mass2_);
    }
    std::string operatorString() const override { return "pocket_potential"; }
};
struct PocketPars : DefaultModelPars {
    static constexpr size_t NScalars=2;
    using FloatType=double;
};
class PocketModel;
using PocketBase=MakeModel(PocketModel,PocketPars);
class PocketModel : public PocketBase {
public:
    double lambda,deltaV,mass2;
    PocketModel(ParameterParser& parser,RunParameters<double>& run,auto toolbox,
                const std::vector<double>& state)
      : PocketBase(parser,run.getLatParams(),toolbox,run.dt,"pockets") {
        alpha=fStar=omegaStar=1.;
        lambda=parser.get<double>("lambda")();
        deltaV=parser.get<double>("delta_v")();
        const double massTrue=parser.get<double>("mass_tv")();
        mass2=massTrue*massTrue;
        if (!(lambda>0.) || !(deltaV>=0.) || !(massTrue>=0.))
            throw std::invalid_argument("Invalid pocket parameters");
        auto phi=fldS(0_c).getFullNDHostView(), ax=fldS(1_c).getFullNDHostView();
        auto pphi=piS(0_c).getFullNDHostView(), pax=piS(1_c).getFullNDHostView();
        const auto layout=fldS(0_c).getLayout();
        for (ptrdiff_t i=0;i<ptrdiff_t(phi.extent(0));++i)
          for (ptrdiff_t j=0;j<ptrdiff_t(phi.extent(1));++j)
            for (ptrdiff_t k=0;k<ptrdiff_t(phi.extent(2));++k) {
                device::IdxArray<3> xyz;
                layout.putSpatialLocationFromMemoryIndexInto(xyz,i,j,k);
                for (auto& v:xyz) v=(v%run.N+run.N)%run.N;
                const auto p=4*((xyz[0]*run.N+xyz[1])*run.N+xyz[2]);
                phi(i,j,k)=state[p];ax(i,j,k)=state[p+1];pphi(i,j,k)=state[p+2];pax(i,j,k)=state[p+3];
            }
        fldS(0_c).getMemoryManager()->pushHostView();fldS(1_c).getMemoryManager()->pushHostView();
        piS(0_c).getMemoryManager()->pushHostView();piS(1_c).getMemoryManager()->pushHostView();
        fldS(0_c).setGhostsAreStale();fldS(1_c).setGhostsAreStale();
        piS(0_c).setGhostsAreStale();piS(1_c).setGhostsAreStale();
        fldS0={1.,0.};piS0={0.,0.};
    }
    auto potentialTerms(Tag<0>) {
        auto phi=fldS(0_c),ax=fldS(1_c);
        return PocketFunction<decltype(phi),decltype(ax),0>(phi,ax,lambda,deltaV,mass2);
    }
    template<int N> auto potDeriv(Tag<N>) {
        auto phi=fldS(0_c),ax=fldS(1_c);
        return PocketFunction<decltype(phi),decltype(ax),N+1>(phi,ax,lambda,deltaV,mass2);
    }
};
}
#endif
