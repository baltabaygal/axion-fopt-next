#ifndef POCKETS_FORMATION_POTENTIAL_H
#define POCKETS_FORMATION_POTENTIAL_H

#include <cmath>
#include <initializer_list>
#include <stdexcept>

// Canonical excursion-unit dilaton: false vacuum phi=0, true vacuum phi=1.
// Axion dependence follows arXiv:2407.01676 PDF Eq. (2.2). U is an explicit
// tilted-quartic choice, not a uniquely specified vacuum model in that paper.
// Optional CPU lattice model selection uses this kernel. Framework evolution
// validation is required before interpreting its runs as scientific evidence.
namespace pockets {
struct FormationParameters {
    double lambda=8., delta_v=.01, mass_false=.1;
    double kappa=4.605170185988092, fa=10.;
    void validate() const {
        if (!std::isfinite(lambda) || !std::isfinite(delta_v) ||
            !std::isfinite(mass_false) || !std::isfinite(kappa) ||
            !std::isfinite(fa) || lambda<=0 || delta_v<0 ||
            mass_false<=0 || kappa<0 || fa<=0)
            throw std::invalid_argument("Invalid formation parameters");
    }
};
struct FormationEvaluation {
    double vacuum, axion, mass_squared;
    double dphi, da, hphi_phi, hphi_a, ha_a;
    double axion_dphi;
    double potential() const { return vacuum+axion; }
    // Positive Q adds to axion-sector energy: dE_a/dt + div(S_a) = Q.
    double axion_exchange(double phi_dot) const { return axion_dphi*phi_dot; }
};
inline FormationEvaluation formationPotential(double phi, double a,
                                             const FormationParameters& p) {
    p.validate();
    if (!std::isfinite(phi) || !std::isfinite(a))
        throw std::invalid_argument("Require finite formation fields");
    const double q=1.-phi, theta=a/p.fa;
    const double m2=p.mass_false*p.mass_false*std::exp(p.kappa*phi);
    const double s=std::sin(theta/2.);
    const double va=2.*p.fa*p.fa*m2*s*s;
    // Factored vacuum energy avoids cancellation close to the true vacuum.
    const double u=q*q*(p.lambda/4.*phi*phi+p.delta_v*(1.+2.*phi+3.*phi*phi));
    const double up=p.lambda/2.*phi*q*(1.-2.*phi)+12.*p.delta_v*phi*phi*(phi-1.);
    const double upp=p.lambda/2.*(1.-6.*phi+6.*phi*phi)+12.*p.delta_v*(3.*phi*phi-2.*phi);
    const double vap=p.kappa*va, aa=p.fa*m2*std::sin(theta);
    FormationEvaluation out{u,va,m2,up+vap,aa,upp+p.kappa*vap,
                            p.kappa*aa,m2*std::cos(theta),vap};
    for (double value : {out.vacuum,out.axion,out.mass_squared,out.dphi,out.da,
                         out.hphi_phi,out.hphi_a,out.ha_a,out.axion_dphi})
        if (!std::isfinite(value))
            throw std::overflow_error("Nonfinite formation potential; no clipping applied");
    return out;
}
} // namespace pockets
#endif
