#ifndef AXION_FOPT_NEXT_ARRIVAL_MAP_HPP
#define AXION_FOPT_NEXT_ARRIVAL_MAP_HPP
// Standalone host-side transfer/force prototype. No TempLat integration yet.
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace axion_fopt_next {
class ArrivalMap {
    std::uint32_t n_;
    double box_;
    std::vector<double> arrival_;

    static std::uint64_t read_le(std::istream& in, unsigned count) {
        std::uint64_t value = 0;
        for (unsigned i = 0; i < count; ++i) {
            const int c = in.get();
            if (c == EOF) throw std::runtime_error("Truncated arrival map");
            value |= std::uint64_t(static_cast<unsigned char>(c)) << (8 * i);
        }
        return value;
    }
    static double read_double(std::istream& in) {
        static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
        return std::bit_cast<double>(read_le(in, 8));
    }
public:
    explicit ArrivalMap(const std::string& path, std::size_t max_cells = 16777216) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open arrival map");
        std::array<char, 8> magic{};
        in.read(magic.data(), magic.size());
        if (std::string(magic.data(), magic.size()) != "AFOPT001")
            throw std::runtime_error("Unsupported arrival-map format");
        n_ = static_cast<std::uint32_t>(read_le(in, 4));
        box_ = read_double(in);
        if (!n_ || !std::isfinite(box_) || box_ <= 0)
            throw std::runtime_error("Invalid arrival-map grid");
        // Divide before multiplying to prevent overflow and unbounded allocation.
        std::size_t cells = n_;
        for (int i = 0; i < 2; ++i) {
            if (cells > max_cells / n_) throw std::runtime_error("Arrival map exceeds cell budget");
            cells *= n_;
        }
        if (cells > max_cells) throw std::runtime_error("Arrival map exceeds cell budget");
        arrival_.reserve(cells);
        for (std::size_t i = 0; i < cells; ++i) {
            double t = read_double(in);
            if (std::isnan(t) || t == -std::numeric_limits<double>::infinity())
                throw std::runtime_error("Invalid arrival time");
            arrival_.push_back(t);
        }
        if (in.get() != EOF) throw std::runtime_error("Trailing arrival-map data");
    }
    std::uint32_t n() const { return n_; }
    double box_size() const { return box_; }
    std::size_t size() const { return arrival_.size(); }
    double arrival(std::size_t flat_index) const { return arrival_.at(flat_index); }
    std::size_t index(std::uint32_t x, std::uint32_t y, std::uint32_t z) const {
        if (x >= n_ || y >= n_ || z >= n_) throw std::out_of_range("Grid index");
        return (std::size_t(x) * n_ + y) * n_ + z;
    }
    double mass_squared(std::size_t cell, double stage_time, double tau_origin,
                        double omega_star, double mass_false, double mass_true) const {
        if (!std::isfinite(stage_time) || !std::isfinite(tau_origin) ||
            !std::isfinite(omega_star) || omega_star <= 0 ||
            !std::isfinite(mass_false) || mass_false < 0 ||
            !std::isfinite(mass_true) || mass_true < 0)
            throw std::invalid_argument("Invalid time or mass units");
        // alpha=1: tau = tau_origin + program_time / omegaStar.
        const double tau = tau_origin + stage_time / omega_star;
        if (!std::isfinite(tau)) throw std::invalid_argument("Conformal time overflow");
        const double mass = tau >= arrival(cell) ? mass_true : mass_false;
        const double scaled = mass / omega_star;
        const double result = scaled * scaled;
        if (!std::isfinite(result)) throw std::invalid_argument("Program mass overflow");
        return result;
    }
};
} // namespace axion_fopt_next
#endif
