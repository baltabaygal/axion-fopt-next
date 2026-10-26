#ifndef STOCHASTIC_V2_MODEL_H
#define STOCHASTIC_V2_MODEL_H
#include "TempLat/util/staticif.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "CosmoInterface/abstractmodel.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include "potential.h"

namespace TempLat {
template<class X, class Y, int Order>
class BubbleFunction : public BinaryOperator<X,Y> {
    double lambda_,dv_,m_,f_;
public:
    BubbleFunction(const X& x,const Y& y,double l,double d,double m,double f)
      : BinaryOperator<X,Y>(x,y),lambda_(l),dv_(d),m_(m),f_(f) {}
    template<typename... IDX>
      requires requires(std::decay_t<X> x,std::decay_t<Y> y,IDX... idx) {
        requires IsVariadicIndex<IDX...>; DoEval::eval(x,idx...); DoEval::eval(y,idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        return bubbleValue<Order>(DoEval::eval(this->mR,idx...),DoEval::eval(this->mT,idx...),lambda_,dv_,m_,f_);
    }
    std::string operatorString() const override { return "pocket_potential"; }
};
struct BubblePars : DefaultModelPars {
    static constexpr size_t NScalars=2;
    using FloatType=double;
};
class BubbleModel;
using BubbleBase=MakeModel(BubbleModel,BubblePars);
class BubbleModel : public BubbleBase {
public:
    double lambda,deltaV,mass,fa;
    BubbleModel(ParameterParser& parser,RunParameters<double>& run,auto toolbox,
                const auto& state)
      : BubbleBase(parser,run.getLatParams(),toolbox,run.dt,"dynamical_bubble") {
        alpha=fStar=omegaStar=1.;
        lambda=parser.get<double>("lambda")(); deltaV=parser.get<double>("delta_v")();
        mass=parser.get<double>("mass_tv")(); fa=parser.get<double>("axion_fa")();
        for(double v:{lambda,deltaV,mass,fa})if(!std::isfinite(v))throw std::invalid_argument("Nonfinite potential parameter");
        if(lambda<=0||deltaV<0||mass<0||fa<=0)throw std::invalid_argument("Invalid potential parameter");
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
        return BubbleFunction<decltype(phi),decltype(ax),0>(phi,ax,lambda,deltaV,mass,fa);
    }
    template<int N> auto potDeriv(Tag<N>) {
        auto phi=fldS(0_c),ax=fldS(1_c);
        return BubbleFunction<decltype(phi),decltype(ax),N+1>(phi,ax,lambda,deltaV,mass,fa);
    }
};
}
#endif
