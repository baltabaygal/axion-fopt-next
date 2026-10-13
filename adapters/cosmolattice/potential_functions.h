#ifndef AXION_FOPT_NEXT_POTENTIAL_FUNCTIONS_H
#define AXION_FOPT_NEXT_POTENTIAL_FUNCTIONS_H
#include "TempLat/lattice/algebra/operators/unaryoperator.h"

namespace TempLat {
template<int Order>
DEVICE_INLINE_FUNCTION double axionPotentialValue(double theta, bool quadratic) {
    if constexpr (Order == 0) return quadratic ? .5 * theta * theta : 1. - device::cos(theta);
    if constexpr (Order == 1) return quadratic ? theta : device::sin(theta);
    if constexpr (Order == 2) return quadratic ? 1. : device::cos(theta);
}

// One expression type for both potentials, without evaluating unused trig functions.
template<typename R, int Order> class AxionPotentialFunction : public UnaryOperator<R> {
    bool quadratic_;
public:
    AxionPotentialFunction(const R& field, bool quadratic) : UnaryOperator<R>(field), quadratic_(quadratic) {}
    template<typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        return axionPotentialValue<Order>(DoEval::eval(this->mR, idx...), quadratic_);
    }
    std::string operatorString() const override { return "axion_potential"; }
};

// At |z|>=20, CPU double-precision tanh already rounds to +/-1.
// Keep the original formula throughout the resolved wall and its tails.
DEVICE_INLINE_FUNCTION double onePlusWallTanhValue(double z) {
    if (z >= 20.) return 2.;
    if (z <= -20.) return 0.;
    return 1. + device::tanh(z);
}
template<typename R> class OnePlusWallTanh : public UnaryOperator<R> {
public:
    explicit OnePlusWallTanh(const R& argument) : UnaryOperator<R>(argument) {}
    template<typename... IDX>
      requires requires(std::decay_t<R> r, IDX... idx) {
        requires IsVariadicIndex<IDX...>;
        DoEval::eval(r, idx...);
      }
    DEVICE_INLINE_FUNCTION double eval(const IDX&... idx) const {
        const double z = DoEval::eval(this->mR, idx...);
        return onePlusWallTanhValue(z);
    }
    std::string operatorString() const override { return "one_plus_wall_tanh"; }
};
} // namespace TempLat
#endif
