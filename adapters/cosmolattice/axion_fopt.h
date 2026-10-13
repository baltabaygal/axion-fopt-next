#ifndef AXION_FOPT_NEXT_COSMOLATTICE_MODEL_H
#define AXION_FOPT_NEXT_COSMOLATTICE_MODEL_H
#include "TempLat/util/staticif.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "CosmoInterface/abstractmodel.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "arrival_map.hpp"
#include "potential_functions.h"
#include <chrono>
#include <unordered_map>

namespace TempLat {
struct AxionFoptPars : DefaultModelPars {
    static constexpr size_t NScalars = 1;
    static constexpr size_t NPotTerms = 1;
    using FloatType = double;
};
class AxionFopt;
using AxionFoptBase = MakeModel(AxionFopt, AxionFoptPars);

class AxionFopt : public AxionFoptBase {
    double tauOrigin_, mf2_, mt2_, wallWidth_, wallSpeed_, quadratic_;
    bool profiling_;
    double massUpdateSeconds_ = 0.;
    double lastTau_ = std::numeric_limits<double>::quiet_NaN(), pendingU_ = 0.;
    bool externalDirty_ = false;
    std::vector<double> arrivalLevels_, massLevels_;
    std::vector<std::uint32_t> massLevelIndices_;
public:
    Field<double, 3> arrival;
    Field<double, 3> mass2;

    AxionFopt(ParameterParser& parser, RunParameters<double>& run,
              auto toolbox, const axion_fopt_next::ArrivalMap& map)
        : AxionFoptBase(parser, run.getLatParams(), toolbox, run.dt, "axion_fopt"),
          tauOrigin_(parser.get<double>("tau_origin", 0.)()),
          arrival("external_arrival", toolbox, run.getLatParams()),
          mass2("external_mass2", toolbox, run.getLatParams()) {
        alpha = 1.;
        fStar = parser.get<double>("f_axion", 1.)();
        omegaStar = parser.get<double>("omega_star", 1.)();
        const double mf = parser.get<double>("mass_false", 0.)();
        const double mt = parser.get<double>("mass_true", 1.)();
        if (!std::isfinite(fStar) || fStar <= 0. || !std::isfinite(omegaStar) || omegaStar <= 0. ||
            !std::isfinite(tauOrigin_) || !std::isfinite(mf) || mf < 0. || !std::isfinite(mt) || mt < 0.)
            throw std::invalid_argument("Invalid axion model units or masses");
        mf2_ = (mf / omegaStar) * (mf / omegaStar);
        mt2_ = (mt / omegaStar) * (mt / omegaStar);
        if (!std::isfinite(mf2_) || !std::isfinite(mt2_)) throw std::invalid_argument("Mass overflow");
        wallWidth_ = parser.get<double>("wall_width", 0.)();
        wallSpeed_ = parser.get<double>("profile_wall_speed", 1.)();
        const auto potential = parser.get<std::string>("potential", "cosine")();
        if (potential != "cosine" && potential != "quadratic")
            throw std::invalid_argument("Require potential=cosine or quadratic");
        quadratic_ = potential == "quadratic" ? 1. : 0.;
        profiling_ = parser.get<bool>("profile", false)();
        if (!std::isfinite(wallWidth_) || wallWidth_ < 0. ||
            !std::isfinite(wallSpeed_) || wallSpeed_ <= 0. || wallSpeed_ > 1.)
            throw std::invalid_argument("Invalid wall width or profile wall speed");
        if (run.N != map.n() || std::abs(run.lSide - omegaStar * map.box_size()) > 1e-12 * run.lSide)
            throw std::invalid_argument("Arrival-map N or program box length mismatch");
        const double theta = parser.get<double>("theta_initial", .4)();
        const double derivative = parser.get<double>("theta_prime_initial", 0.)();
        if (!std::isfinite(theta) || !std::isfinite(derivative))
            throw std::invalid_argument("Invalid homogeneous initial conditions");
        // Initialize explicitly, without vacuum fluctuations or Friedmann feedback.
        fldS(0_c) = theta;
        piS(0_c) = derivative / omegaStar; // a_initial=1; pi=a^2*dtheta/du.
        fldS0 = {fStar * theta};
        piS0 = {fStar * derivative};

        auto host = arrival.getFullNDHostView();
        const auto layout = arrival.getLayout();
        bool cache = wallWidth_ > 0. && parser.get<bool>("mass_cache", true)();
        std::unordered_map<double, std::uint32_t> levels;
        const std::size_t maxLevels = std::max<std::size_t>(1, map.size() / 8);
        if (cache) massLevelIndices_.reserve(host.extent(0)*host.extent(1)*host.extent(2));
        for (ptrdiff_t i = 0; i < ptrdiff_t(host.extent(0)); ++i)
          for (ptrdiff_t j = 0; j < ptrdiff_t(host.extent(1)); ++j)
            for (ptrdiff_t k = 0; k < ptrdiff_t(host.extent(2)); ++k) {
                device::IdxArray<3> xyz;
                layout.putSpatialLocationFromMemoryIndexInto(xyz, i, j, k);
                for (auto& x : xyz) x = (x % run.N + run.N) % run.N;
                host(i, j, k) = map.arrival(map.index(xyz[0], xyz[1], xyz[2]));
                if (cache) {
                    const double value = host(i,j,k);
                    const auto [level, inserted] = levels.try_emplace(value, static_cast<std::uint32_t>(arrivalLevels_.size()));
                    if (inserted) arrivalLevels_.push_back(value);
                    massLevelIndices_.push_back(level->second);
                    if (arrivalLevels_.size() > maxLevels) {
                        // Unstructured maps fall back to the ordinary field expression.
                        cache = false;
                        levels.clear();
                        arrivalLevels_.clear(); arrivalLevels_.shrink_to_fit();
                        massLevelIndices_.clear(); massLevelIndices_.shrink_to_fit();
                    }
                }
            }
        massLevels_.resize(arrivalLevels_.size());
        arrival.getMemoryManager()->pushHostView();
        arrival.setGhostsAreStale();
        setExternalTime(0.);
    }

    double conformalTime(double u) const { return tauOrigin_ + u / omegaStar; }
    double massUpdateSeconds() const { return massUpdateSeconds_; }
    void deferExternalTime(double u) {
        if (!std::isfinite(conformalTime(u))) throw std::invalid_argument("Nonfinite conformal endpoint time");
        pendingU_ = u;
        externalDirty_ = true;
    }
    void ensureExternalMass() { if (externalDirty_) setExternalTime(pendingU_); }
    void setExternalTime(double u) {
        const double tau = conformalTime(u);
        if (!std::isfinite(tau)) throw std::invalid_argument("Nonfinite conformal stage time");
        externalDirty_ = false;
        if (tau == lastTau_) return;
        const auto begin = profiling_ ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        if (!massLevelIndices_.empty()) {
            const double scale = wallSpeed_ / wallWidth_;
            for (std::size_t i = 0; i < arrivalLevels_.size(); ++i)
                massLevels_[i] = mf2_ + (mt2_ - mf2_) * .5 * onePlusWallTanhValue((tau-arrivalLevels_[i])*scale);
            auto host = mass2.getFullNDHostView();
            std::size_t index = 0;
            for (ptrdiff_t i = 0; i < ptrdiff_t(host.extent(0)); ++i)
              for (ptrdiff_t j = 0; j < ptrdiff_t(host.extent(1)); ++j)
                for (ptrdiff_t k = 0; k < ptrdiff_t(host.extent(2)); ++k)
                    host(i,j,k) = massLevels_[massLevelIndices_[index++]];
            mass2.getMemoryManager()->pushHostView();
            mass2.setGhostsAreStale();
        } else if (wallWidth_ == 0.)
            mass2 = mf2_ + (mt2_ - mf2_) * heaviside(tau - arrival);
        else {
            auto argument = (tau - arrival) * (wallSpeed_ / wallWidth_);
            mass2 = mf2_ + (mt2_ - mf2_) * .5 *
                OnePlusWallTanh<decltype(argument)>(argument);
        }
        lastTau_ = tau;
        if (profiling_) massUpdateSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    }
    auto potentialShape() const {
        auto field = fldS(0_c);
        return AxionPotentialFunction<decltype(field), 0>(field, quadratic_ != 0.);
    }
    auto potentialTerms(Tag<0>) { ensureExternalMass(); return mass2 * potentialShape(); }
    auto potDeriv(Tag<0>) {
        ensureExternalMass();
        auto field = fldS(0_c);
        return mass2 * AxionPotentialFunction<decltype(field), 1>(field, quadratic_ != 0.);
    }
    auto potDeriv2(Tag<0>) {
        ensureExternalMass();
        auto field = fldS(0_c);
        return mass2 * AxionPotentialFunction<decltype(field), 2>(field, quadratic_ != 0.);
    }
    double shapeValue(double theta) const { return axionPotentialValue<0>(theta, quadratic_ != 0.); }
    double externalPowerRate() const {
        return wallWidth_ > 0. && mt2_ != mf2_
            ? 2. * wallSpeed_ / (wallWidth_ * omegaStar * (mt2_ - mf2_)) : 0.;
    }
    double falseMassSquared() const { return mf2_; }
    double trueMassSquared() const { return mt2_; }
    auto externalPower() {
        ensureExternalMass();
        // d(m_program^2)/du times the potential shape. Valid also for mf=mt.
        const double rate = externalPowerRate();
        return rate * (mass2 - mf2_) * (mt2_ - mass2) * potentialShape();
    }
};
} // namespace TempLat
#endif
