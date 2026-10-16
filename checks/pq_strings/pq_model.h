#ifndef PQ_STRINGS_MODEL_H
#define PQ_STRINGS_MODEL_H
#include "TempLat/util/staticif.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "CosmoInterface/abstractmodel.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include <cmath>
#include <vector>

namespace TempLat {
// f_a=omegaStar=fStar=1. Order 0: V; 1,2: Cartesian derivatives.
template<int Order>
DEVICE_INLINE_FUNCTION double pqValue(double x, double y, double lambda, double chi) {
    const double r2=x*x+y*y, r=device::sqrt(r2);
    if constexpr (Order==0) return .25*lambda*(r2-1.)*(r2-1.)+chi*(r2-x*r);
    if (r==0.) return 0.; // Unique continuous gradient limit; no epsilon.
    if constexpr (Order==1) return lambda*(r2-1.)*x+chi*(2.*x-r-x*x/r);
    if constexpr (Order==2) return lambda*(r2-1.)*y+chi*(2.*y-x*y/r);
}
template<class X, class Y, int Order>
class PQFunction : public BinaryOperator<X,Y> {
    double lambda_,chi_;
public:
    PQFunction(const X& x,const Y& y,double lambda,double chi)
      : BinaryOperator<X,Y>(x,y),lambda_(lambda),chi_(chi) {}
    template<typename... IDX>
      requires requires(std::decay_t<X> x,std::decay_t<Y> y,IDX... idx) {
        requires IsVariadicIndex<IDX...>; DoEval::eval(x,idx...); DoEval::eval(y,idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        return pqValue<Order>(DoEval::eval(this->mR,idx...),DoEval::eval(this->mT,idx...),lambda_,chi_);
    }
    std::string operatorString() const override { return "pq_potential"; }
};
struct PQPars : DefaultModelPars {
    static constexpr size_t NScalars=2;
    using FloatType=double;
};
class PQModel;
using PQBase=MakeModel(PQModel,PQPars);
class PQModel : public PQBase {
public:
    double lambda,lambdaPhysical,lambdaDerivative=0.,coreGrowthUntil=0.;
    double tauOrigin,tauStar,exponent,chiStar,chi=0.,chiDerivative=0.;
    std::string massSchedule;
    double massOnset,massRamp,massCap;
    PQModel(ParameterParser& parser,RunParameters<double>& run,auto toolbox,
            const std::vector<double>& state)
      : PQBase(parser,run.getLatParams(),toolbox,run.dt,"pq_strings") {
        alpha=fStar=omegaStar=1.;
        lambdaPhysical=parser.get<double>("lambda")();
        coreGrowthUntil=parser.get<double>("core_growth_until",0.)();
        tauOrigin=parser.get<double>("tau_origin")();
        tauStar=parser.get<double>("tau_star")();
        exponent=parser.get<double>("n_tau",8.)();
        chiStar=parser.get<double>("chi_star",std::pow(tauOrigin/(tauStar*tauStar),2))();
        massSchedule=parser.get<std::string>("mass_schedule","powerlaw")();
        massOnset=parser.get<double>("mass_onset",tauOrigin)();
        massRamp=parser.get<double>("mass_ramp_duration",1.)();
        massCap=parser.get<double>("mass_cap",0.)();
        if ((massSchedule!="powerlaw" && massSchedule!="smooth_cap") ||
            massOnset<tauOrigin || !(massRamp>0.) || massCap<0.)
            throw std::invalid_argument("Invalid susceptibility schedule");
        if (!(lambdaPhysical>0.) || !(tauOrigin>0.) || !(tauStar>0.) || exponent<0. || chiStar<0. ||
            (coreGrowthUntil!=0. && coreGrowthUntil<=tauOrigin))
            throw std::invalid_argument("Invalid PQ parameters");
        auto x=fldS(0_c).getFullNDHostView(), y=fldS(1_c).getFullNDHostView();
        auto px=piS(0_c).getFullNDHostView(), py=piS(1_c).getFullNDHostView();
        const auto layout=fldS(0_c).getLayout();
        for (ptrdiff_t i=0;i<ptrdiff_t(x.extent(0));++i)
          for (ptrdiff_t j=0;j<ptrdiff_t(x.extent(1));++j)
            for (ptrdiff_t k=0;k<ptrdiff_t(x.extent(2));++k) {
                device::IdxArray<3> xyz;
                layout.putSpatialLocationFromMemoryIndexInto(xyz,i,j,k);
                for (auto& v:xyz) v=(v%run.N+run.N)%run.N;
                const auto p=4*((xyz[0]*run.N+xyz[1])*run.N+xyz[2]);
                x(i,j,k)=state[p];y(i,j,k)=state[p+1];px(i,j,k)=state[p+2];py(i,j,k)=state[p+3];
            }
        fldS(0_c).getMemoryManager()->pushHostView();fldS(1_c).getMemoryManager()->pushHostView();
        piS(0_c).getMemoryManager()->pushHostView();piS(1_c).getMemoryManager()->pushHostView();
        fldS(0_c).setGhostsAreStale();fldS(1_c).setGhostsAreStale();
        piS(0_c).setGhostsAreStale();piS(1_c).setGhostsAreStale();
        fldS0={1.,0.};piS0={0.,0.};
        setExternalTime(0.);
    }
    double conformalTime(double u) const { return tauOrigin+u; }
    void setExternalTime(double u) {
        const double tau=conformalTime(u);
        if (!(tau>0.) || !std::isfinite(tau)) throw std::invalid_argument("Invalid stage time");
        // Preparation only: s=-1 gives lambda proportional to a^-4 and w_core proportional to tau.
        // Match continuously onto constant physical lambda at tau_cg.
        const bool growing=coreGrowthUntil>0. && tau<coreGrowthUntil;
        lambda=growing ? lambdaPhysical*std::pow(coreGrowthUntil/tau,4) : lambdaPhysical;
        lambdaDerivative=growing ? -4.*lambda/tau : 0.;
        if (massSchedule=="smooth_cap") {
            const double s=std::clamp((tau-massOnset)/massRamp,0.,1.);
            chi=massCap*massCap*s*s*s*(10.+s*(-15.+6.*s));
            chiDerivative=massCap*massCap*30.*s*s*(s-1.)*(s-1.)/massRamp;
        } else {
            chi=chiStar*std::pow(tau/tauStar,exponent);
            chiDerivative=exponent*chi/tau;
        }
        if (!std::isfinite(chi) || !std::isfinite(chiDerivative) || !std::isfinite(lambda) || !std::isfinite(lambdaDerivative))
            throw std::invalid_argument("External coefficient overflow");
    }
    auto potentialTerms(Tag<0>) {
        auto x=fldS(0_c),y=fldS(1_c);
        return PQFunction<decltype(x),decltype(y),0>(x,y,lambda,chi);
    }
    template<int N> auto potDeriv(Tag<N>) {
        auto x=fldS(0_c),y=fldS(1_c);
        return PQFunction<decltype(x),decltype(y),N+1>(x,y,lambda,chi);
    }
};
}
#endif
