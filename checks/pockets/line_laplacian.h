#ifndef POCKETS_LINE_LAPLACIAN_H
#define POCKETS_LINE_LAPLACIAN_H
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
namespace TempLat {
// Canonical potential derivatives stay untouched. This hook is used only
// by the explicitly generated, user-owned scalar-kernel overlay.
template<class R> class PocketLineLaplacian : public UnaryOperator<R> {
    using UnaryOperator<R>::mR;
    double dx2_;int order_;
public:
    PocketLineLaplacian(R field,int order):UnaryOperator<R>(field),dx2_(std::pow(GetDx::getDx(field),2)),order_(order) {
        static_assert(GetNDim::get<R>()==1,"Planar Laplacian requires NDim=1");
    }
    void doWeNeedGhosts() const {mR.confirmGhostsUpToDate();}
    template<class I> requires IsVariadicIndex<I>
    DEVICE_INLINE_FUNCTION auto eval(const I& i) const {
        const auto centre=DoEval::eval(mR,i);
        auto near=-2*centre;near+=DoEval::eval(mR,i+1);near+=DoEval::eval(mR,i-1);
        if(order_==2)return near/dx2_;
        auto far=-2*centre;far+=DoEval::eval(mR,i+2);far+=DoEval::eval(mR,i-2);
        return ((4./3.)*near-(1./12.)*far)/dx2_;
    }
    std::string operatorString() const override {return "planar_laplacian";}
};
template<class Model,int N> auto PocketPlanarSpatialLaplacian(Model& model,Tag<N> n) {
    if constexpr(requires {model.lineSpatialLaplacian(n);})return model.lineSpatialLaplacian(n);
    else return LatLapl(model.fldS(n));
}
}
#endif
