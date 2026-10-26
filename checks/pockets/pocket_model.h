#ifndef POCKETS_MODEL_H
#define POCKETS_MODEL_H
#include "TempLat/util/staticif.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "CosmoInterface/abstractmodel.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include "pocket_potential.h"

namespace TempLat {
template<class X, class Y, int Order>
class PocketFunction : public BinaryOperator<X,Y> {
    double lambda_,deltaV_,mass2_,stiffening_;
    int power_;
    bool useFormation_;
    pockets::FormationParameters formation_;
public:
    PocketFunction(const X& x,const Y& y,double lambda,double deltaV,double mass2,int power,double stiffening,const pockets::FormationParameters* formation=nullptr)
      : BinaryOperator<X,Y>(x,y),lambda_(lambda),deltaV_(deltaV),mass2_(mass2),stiffening_(stiffening),power_(power),useFormation_(formation!=nullptr),formation_(formation?*formation:pockets::FormationParameters{}) {}
    template<typename... IDX>
      requires requires(std::decay_t<X> x,std::decay_t<Y> y,IDX... idx) {
        requires IsVariadicIndex<IDX...>; DoEval::eval(x,idx...); DoEval::eval(y,idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        return pocketValue<Order>(DoEval::eval(this->mR,idx...),DoEval::eval(this->mT,idx...),lambda_,deltaV_,mass2_,power_,stiffening_,useFormation_?&formation_:nullptr);
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
    double lambda,deltaV,mass2,vacuumStiffening;
    int massPower;
    bool useFormation;
    pockets::FormationParameters formation;
    const pockets::FormationParameters* formationParameters() const {return useFormation?&formation:nullptr;}
    PocketModel(ParameterParser& parser,RunParameters<double>& run,auto toolbox,
                const auto& state)
      : PocketBase(parser,run.getLatParams(),toolbox,run.dt,"pockets") {
        alpha=fStar=omegaStar=1.;
        lambda=parser.get<double>("lambda")();
        deltaV=parser.get<double>("delta_v")();
        const double massTrue=parser.get<double>("mass_tv")();
        mass2=massTrue*massTrue;
        massPower=parser.get<int>("mass_power",2)();
        vacuumStiffening=parser.get<double>("vacuum_stiffening",0.)();
        if(!(vacuumStiffening>=0. && vacuumStiffening<16.))throw std::invalid_argument("Require 0<=vacuum_stiffening<16");
        if(massPower!=2 && massPower!=4) throw std::invalid_argument("mass_power must be 2 or 4");
        if (!(lambda>0.) || !(deltaV>=0.) || !(massTrue>=0.))
            throw std::invalid_argument("Invalid pocket parameters");
        const auto potentialModel=parser.get<std::string>("potential_model","polynomial")();
        if(potentialModel!="polynomial" && potentialModel!="formation")
            throw std::invalid_argument("potential_model must be polynomial or formation");
        useFormation=potentialModel=="formation";
        if(useFormation) {
            if(vacuumStiffening!=0. || massPower!=2)
                throw std::invalid_argument("Formation model cannot use polynomial stiffening or mass_power");
            formation={lambda,deltaV,parser.get<double>("mass_fv")(),
                       parser.get<double>("mass_kappa")(),parser.get<double>("axion_fa")()};
            formation.validate();
            const double predicted=formation.mass_false*std::exp(formation.kappa/2.);
            if(!std::isfinite(predicted) || std::abs(predicted-massTrue)>1e-10*std::max(predicted,massTrue))
                throw std::invalid_argument("mass_tv must equal mass_fv*exp(mass_kappa/2)");
        }
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
        return PocketFunction<decltype(phi),decltype(ax),0>(phi,ax,lambda,deltaV,mass2,massPower,vacuumStiffening,formationParameters());
    }
    template<int N> auto potDeriv(Tag<N>) {
        auto phi=fldS(0_c),ax=fldS(1_c);
        return PocketFunction<decltype(phi),decltype(ax),N+1>(phi,ax,lambda,deltaV,mass2,massPower,vacuumStiffening,formationParameters());
    }
};
}
#endif
