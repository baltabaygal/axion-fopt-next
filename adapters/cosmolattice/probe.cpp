#include "arrival_map.hpp"
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc != 8) throw std::invalid_argument("Usage: probe map stage_time tau_origin omega_star theta mass_false mass_true");
        axion_fopt_next::ArrivalMap map(argv[1]);
        const double stage = std::stod(argv[2]), origin = std::stod(argv[3]);
        const double omega = std::stod(argv[4]), theta = std::stod(argv[5]);
        const double mf = std::stod(argv[6]), mt = std::stod(argv[7]);
        if (!std::isfinite(theta)) throw std::invalid_argument("Invalid angle");
        std::cout << std::setprecision(17);
        for (std::size_t cell = 0; cell < map.size(); ++cell) {
            const double m2 = map.mass_squared(cell, stage, origin, omega, mf, mt);
            std::cout << m2 << ' ' << m2 * std::sin(theta) << ' '
                      << m2 * (1 - std::cos(theta)) << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
