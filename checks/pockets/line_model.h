#ifndef POCKETS_LINE_MODEL_H
#define POCKETS_LINE_MODEL_H
// Reuse the production expression wrapper and scalar potential; no new forces.
#include "pocket_model.h"
#include "line_laplacian.h"
namespace TempLat {
struct PocketLinePars : DefaultModelPars {
    static constexpr size_t NScalars=2;
    static constexpr size_t NDim=1;
    using FloatType=double;
};
class PocketLineModel;
using PocketLineBase=MakeModel(PocketLineModel,PocketLinePars);
class PocketLineModel : public PocketLineBase {
public:
    pockets::FormationParameters formation;
    int spatialOrder;
    PocketLineModel(ParameterParser& parser,RunParameters<double>& run,auto toolbox,
                    const std::vector<double>& state)
      : PocketLineBase(parser,run.getLatParams(),toolbox,run.dt,"pocket_line") {
        alpha=fStar=omegaStar=1.;
        spatialOrder=parser.get<int>("line_spatial_order",2)();
        if(spatialOrder!=2&&spatialOrder!=4)throw std::invalid_argument("line_spatial_order must be 2 or 4");
#ifndef POCKETS_LINE_CUSTOM_LAPLACIAN
        if(spatialOrder!=2)throw std::invalid_argument("This line build has no custom Laplacian hook");
#endif
        formation={parser.get<double>("lambda")(),parser.get<double>("delta_v")(),
                   parser.get<double>("mass_fv")(),parser.get<double>("mass_kappa")(),
                   parser.get<double>("axion_fa")()};
        formation.validate();
        auto phi=fldS(0_c).getFullNDHostView(),a=fldS(1_c).getFullNDHostView();
        auto pp=piS(0_c).getFullNDHostView(),pa=piS(1_c).getFullNDHostView();
        const auto layout=fldS(0_c).getLayout();
        for(ptrdiff_t i=0;i<ptrdiff_t(phi.extent(0));++i) {
            device::IdxArray<1> p;layout.putSpatialLocationFromMemoryIndexInto(p,i);
            const auto j=4*((p[0]%run.N+run.N)%run.N);
            phi(i)=state.at(j);a(i)=state.at(j+1);pp(i)=state.at(j+2);pa(i)=state.at(j+3);
        }
        fldS(0_c).getMemoryManager()->pushHostView();fldS(1_c).getMemoryManager()->pushHostView();
        piS(0_c).getMemoryManager()->pushHostView();piS(1_c).getMemoryManager()->pushHostView();
        fldS(0_c).setGhostsAreStale();fldS(1_c).setGhostsAreStale();
        piS(0_c).setGhostsAreStale();piS(1_c).setGhostsAreStale();
        fldS0={1.,0.};piS0={0.,0.};
    }
    template<int N> auto lineSpatialLaplacian(Tag<N> n) {
        auto field=fldS(n);return PocketLineLaplacian<decltype(field)>(field,spatialOrder);
    }
    auto potentialTerms(Tag<0>) {
        auto phi=fldS(0_c),a=fldS(1_c);
        return PocketFunction<decltype(phi),decltype(a),0>(phi,a,formation.lambda,formation.delta_v,0.,2,0.,&formation);
    }
    template<int N> auto potDeriv(Tag<N>) {
        auto phi=fldS(0_c),a=fldS(1_c);
        return PocketFunction<decltype(phi),decltype(a),N+1>(phi,a,formation.lambda,formation.delta_v,0.,2,0.,&formation);
    }
};
}
#endif
