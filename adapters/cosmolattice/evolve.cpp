#include "TempLat/session/sessionguard.h"
#include "axion_fopt.h"
#include "CosmoInterface/evolvers/positionverlet.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/lattice/algebra/spatialderivatives/normgradientsquare.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <chrono>

using namespace TempLat;

static void validateNumericInputs(const ParameterParser& parser) {
    // Upstream stream conversion can turn malformed/nonfinite text into zero.
    // Validate supplied single-value adapter inputs before typed parsing.
    for (const char* name : {"N", "lSide", "kIR", "dx", "dt", "tMax", "t0", "a0", "H0",
                             "omegaEoS", "f_axion", "omega_star", "mass_false", "mass_true",
                             "tau_origin", "theta_initial", "theta_prime_initial", "output_every",
                             "wall_width", "profile_wall_speed", "snapshot_at", "snapshot_every"}) {
        const auto found = parser.getParams().find(name);
        if (found == parser.getParams().end()) continue;
        std::istringstream input(found->second);
        double value;
        if (!(input >> value) || !std::isfinite(value))
            throw std::invalid_argument(std::string("Invalid numeric input: ") + name);
        input >> std::ws;
        if (!input.eof()) throw std::invalid_argument(std::string("Trailing numeric input: ") + name);
        if ((std::string(name) == "N" || std::string(name) == "output_every" || std::string(name) == "snapshot_every") &&
            (value != std::floor(value) || value > 10000000. || value < (std::string(name) == "snapshot_every" ? 0. : 1.)))
            throw std::invalid_argument(std::string("Invalid integer input: ") + name);
    }
}

static void snapshot(AxionFopt& model, const std::filesystem::path& path, ptrdiff_t n,
                     bool binary = false, double box = 0., double u = 0.) {
    model.ensureExternalMass();
    auto theta = model.fldS(0_c).getFullNDHostView();
    auto momentum = model.piS(0_c).getFullNDHostView();
    auto arrival = model.arrival.getFullNDHostView();
    auto mass = model.mass2.getFullNDHostView();
    const auto layout = model.fldS(0_c).getLayout();
    std::ofstream out(path, binary ? std::ios::binary : std::ios::out);
    std::vector<double> values;
    if (binary) {
        static_assert(std::endian::native == std::endian::little);
        const std::uint32_t dimension = static_cast<std::uint32_t>(n);
        out.write("AFST0001", 8);
        out.write(reinterpret_cast<const char*>(&dimension), 4);
        for (const double value : {box, u, model.aI})
            out.write(reinterpret_cast<const char*>(&value), 8);
        values.reserve(4 * std::size_t(n) * n * n);
    } else out << std::setprecision(17) << "x,y,z,theta,pi,arrival,mass2\n";
    for (ptrdiff_t x = 0; x < n; ++x)
      for (ptrdiff_t y = 0; y < n; ++y)
        for (ptrdiff_t z = 0; z < n; ++z) {
            device::IdxArray<3> idx;
            // v2.0.0's returned ownership flag omits the ghost offset.
            // Check the actual interior memory bounds instead.
            layout.putMemoryIndexFromSpatialLocationInto(idx, x, y, z);
            for (size_t d = 0; d < 3; ++d)
                if (idx[d] < layout.getNGhosts() ||
                    idx[d] >= layout.getNGhosts() + layout.getLocalSizes()[d])
                    throw std::runtime_error("Snapshot index outside owned interior");
            auto i = idx[0], j = idx[1], k = idx[2];
            if (binary) {
                for (const double value : {theta(i,j,k), momentum(i,j,k), arrival(i,j,k), mass(i,j,k)})
                    values.push_back(value);
            } else out << x << ',' << y << ',' << z << ',' << theta(i,j,k) << ',' << momentum(i,j,k)
                       << ',' << arrival(i,j,k) << ',' << mass(i,j,k) << '\n';
        }
    if (binary) out.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double));
    if (!out) throw std::runtime_error("Snapshot write failed");
}

struct Diagnostics {
    double theta = 0., pi = 0., theta2 = 0., potential = 0.;
    double kinetic = 0., gradient = 0., inside = 0., source = 0.;
};

static Diagnostics measure(AxionFopt& model, bool energies, double u) {
    model.ensureExternalMass();
    if (energies) model.fldS(0_c).confirmGhostsUpToDate();
    auto theta = model.fldS(0_c).getFullNDHostView();
    auto momentum = model.piS(0_c).getFullNDHostView();
    auto arrival = model.arrival.getFullNDHostView();
    auto mass = model.mass2.getFullNDHostView();
    const auto layout = model.fldS(0_c).getLayout();
    const auto start = layout.getNGhosts();
    const auto sizes = layout.getLocalSizes();
    const double dx2 = model.fldS(0_c).getDx() * model.fldS(0_c).getDx();
    const double tau = model.conformalTime(u), rate = model.externalPowerRate();
    const double mf = model.falseMassSquared(), mt = model.trueMassSquared();
    const double a2 = model.aI * model.aI;
    const double a6 = a2 * a2 * a2;
    Diagnostics d;
    for (ptrdiff_t i = start; i < start + sizes[0]; ++i)
      for (ptrdiff_t j = start; j < start + sizes[1]; ++j)
        for (ptrdiff_t k = start; k < start + sizes[2]; ++k) {
            const double th = theta(i,j,k), pi = momentum(i,j,k), m = mass(i,j,k);
            const double shape = model.shapeValue(th), potential = m * shape;
            d.theta += th; d.pi += pi; d.theta2 += th * th; d.potential += potential;
            if (energies) {
                const double gx = theta(i+1,j,k)-th, gy = theta(i,j+1,k)-th, gz = theta(i,j,k+1)-th;
                const double gradient = .5 * (gx*gx+gy*gy+gz*gz) / (dx2*a2), kinetic = .5*pi*pi/a6;
                d.kinetic += kinetic; d.gradient += gradient;
                if (tau >= arrival(i,j,k)) d.inside += kinetic+gradient+potential;
                d.source += rate * (m-mf) * (mt-m) * shape;
            }
        }
    const double count = double(sizes[0]) * sizes[1] * sizes[2];
    d.theta /= count; d.pi /= count; d.theta2 /= count; d.potential /= count;
    d.kinetic /= count; d.gradient /= count; d.inside /= count; d.source /= count;
    return d;
}

int main(int argc, char** argv) {
    try {
        using Clock = std::chrono::steady_clock;
        const auto start = Clock::now();
        SessionGuard guard(argc, argv, false);
        ParameterParser parser(argc, argv);
        const bool profile = parser.get<bool>("profile", false)();
        double evolutionSeconds = 0., diagnosticsSeconds = 0., snapshotSeconds = 0.;
        validateNumericInputs(parser);
        RunParameters<double> run(parser);
        if (run.eType != PV2 || run.t0 != 0. || run.a0 != 1. || run.withGWs || run.doResolutionPreserving)
            throw std::invalid_argument("Require PV2, t0=0, a0=1, no GWs or resolution preserving");
        if (run.expansion && !run.fixedBackground)
            throw std::invalid_argument("Prescribed mass requires fixedBackground=true or expansion=false");
        if (run.expansion && (!std::isfinite(run.H0) || run.H0 <= 0. ||
                             !std::isfinite(run.omegaEoS) || run.omegaEoS <= -1./3.))
            throw std::invalid_argument("Require positive H0 and omegaEoS > -1/3");
        if (run.N < 4 || !std::isfinite(run.lSide) || run.lSide <= 0. ||
            !std::isfinite(run.dt) || run.dt <= 0. || run.dt / run.dx > .25 ||
            !std::isfinite(run.tMax) || run.tMax <= 0.)
            throw std::invalid_argument("Invalid grid/time or dt/dx > 0.25");
        const double stepsFloat = run.tMax / run.dt;
        if (!std::isfinite(stepsFloat) || stepsFloat > 10000000. ||
            std::abs(stepsFloat - std::round(stepsFloat)) > 1e-9)
            throw std::invalid_argument("tMax must be an integral number of dt steps, at most 10 million");
        const int steps = static_cast<int>(std::round(stepsFloat));
        const int outputEvery = parser.get<int>("output_every", 1)();
        if (outputEvery <= 0) throw std::invalid_argument("output_every must be positive");
        const auto snapshotFormat = parser.get<std::string>("snapshot_format", "csv")();
        if (snapshotFormat != "csv" && snapshotFormat != "binary")
            throw std::invalid_argument("Require snapshot_format=csv or binary");
        const bool binary = snapshotFormat == "binary";
        const int snapshotEvery = parser.get<int>("snapshot_every", 0)();
        if (snapshotEvery < 0) throw std::invalid_argument("snapshot_every must be nonnegative");
        const double snapshotAt = parser.get<double>("snapshot_at", -1.)();
        if (snapshotAt != -1. && (snapshotAt < 0. || snapshotAt > run.tMax))
            throw std::invalid_argument("snapshot_at must lie in [0,tMax]");
        bool checkpointWritten = false;
        const bool energies = parser.get<bool>("energy_diagnostics", false)();
        const std::filesystem::path output(parser.get<std::string>("result_dir")());
        axion_fopt_next::ArrivalMap map(parser.get<std::string>("arrival_map")());
        auto toolbox = MemoryToolBox<3>::makeShared(run.N, 1);
        AxionFopt model(parser, run, toolbox, map);
        if (std::filesystem::exists(output)) throw std::invalid_argument("result_dir already exists; choose a fresh directory");
        std::filesystem::create_directories(output);
        PositionVerlet<double> integrator(model, run);
        std::ofstream history(output / "history.csv");
        history << std::setprecision(17) << "u,tau,a,theta_mean,pi_mean,theta2_mean,potential_mean";
        if (energies) history << ",kinetic_mean,gradient_mean,total_mean,inside_mean,source_mean";
        history << '\n';
        const double initializationSeconds = std::chrono::duration<double>(Clock::now()-start).count();
        auto snapshotStart = profile ? Clock::now() : Clock::time_point{};
        snapshot(model, output / (binary ? "initial.bin" : "initial.csv"), run.N, binary, run.lSide, 0.);
        if (profile) snapshotSeconds += std::chrono::duration<double>(Clock::now()-snapshotStart).count();
        for (int step = 0; step <= steps; ++step) {
            const double u = step * run.dt;
            if (step % outputEvery == 0 || step == steps) {
                const auto diagnosticsStart = profile ? Clock::now() : Clock::time_point{};
                integrator.sync(model, u);
                const auto d = measure(model, energies, u);
                history << u << ',' << model.conformalTime(u) << ',' << model.aI << ','
                        << d.theta << ',' << d.pi << ',' << d.theta2 << ',' << d.potential;
                if (energies) {
                    history << ',' << d.kinetic << ',' << d.gradient << ',' << d.kinetic + d.gradient + d.potential
                            << ',' << d.inside << ',' << d.source;
                }
                history << '\n';
                if (profile) diagnosticsSeconds += std::chrono::duration<double>(Clock::now()-diagnosticsStart).count();
            }
            if (!checkpointWritten && snapshotAt >= 0. && u >= snapshotAt) {
                snapshotStart = profile ? Clock::now() : Clock::time_point{};
                integrator.sync(model, u);
                snapshot(model, output / (binary ? "checkpoint.bin" : "checkpoint.csv"), run.N, binary, run.lSide, u);
                checkpointWritten = true;
                if (profile) snapshotSeconds += std::chrono::duration<double>(Clock::now()-snapshotStart).count();
            }
            if (snapshotEvery > 0 && step % snapshotEvery == 0) {
                snapshotStart = profile ? Clock::now() : Clock::time_point{};
                integrator.sync(model, u);
                std::ostringstream frameName;
                frameName << "frame_" << std::setw(7) << std::setfill('0') << step << (binary ? ".bin" : ".csv");
                snapshot(model, output / frameName.str(), run.N, binary, run.lSide, u);
                if (profile) snapshotSeconds += std::chrono::duration<double>(Clock::now()-snapshotStart).count();
            }
            if (step < steps) {
                const auto evolutionStart = profile ? Clock::now() : Clock::time_point{};
                integrator.evolve(model, u, false);
                if (profile) evolutionSeconds += std::chrono::duration<double>(Clock::now()-evolutionStart).count();
            }
        }
        snapshotStart = profile ? Clock::now() : Clock::time_point{};
        snapshot(model, output / (binary ? "final.bin" : "final.csv"), run.N, binary, run.lSide, run.tMax);
        if (profile) snapshotSeconds += std::chrono::duration<double>(Clock::now()-snapshotStart).count();
        if (!history) throw std::runtime_error("History write failed");
        std::ofstream parameters(output / "parameters.txt");
        parameters << parser;
        if (!parameters) throw std::runtime_error("Parameter write failed");
        if (profile) {
            std::ofstream timings(output / "timings.json");
            timings << std::setprecision(17) << "{\n"
                    << "  \"initialization_seconds\": " << initializationSeconds << ",\n"
                    << "  \"evolution_seconds\": " << evolutionSeconds << ",\n"
                    << "  \"diagnostics_seconds\": " << diagnosticsSeconds << ",\n"
                    << "  \"snapshots_seconds\": " << snapshotSeconds << ",\n"
                    << "  \"mass_update_seconds\": " << model.massUpdateSeconds() << ",\n"
                    << "  \"total_seconds\": " << std::chrono::duration<double>(Clock::now()-start).count() << "\n}\n";
            if (!timings) throw std::runtime_error("Timing write failed");
        }
        std::cout << "Completed " << steps << " CosmoLattice PV2 steps: " << output << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
