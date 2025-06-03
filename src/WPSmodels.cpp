#include <cmath>
#include "real_type.h"


Real calc_S_qq_Amiet_rozenberg(const Real Ux,const Real omega,const Real rho,const Real tau_w,
    const Real delta,const Real delta_star,const Real theta,const Real dpdx) {
    
    
    const Real nu = 1.51e-5;

    // Eq. 2
    Real omega_bar = (omega * delta_star) / Ux;

    Real D = delta / delta_star;

    // Table 1: beta_c clipped between 0.01 and 30.0
    Real beta_c = std::max(std::min((theta / tau_w) * dpdx, 30.0), 0.01);

    // Cf = tau_w / (0.5 * rho * Ux^2)
    Real C_f = tau_w / (0.5 * rho * Ux * Ux);

    // Friction velocity
    Real u_tau = std::sqrt(tau_w / rho);

    // Eq. 4
    Real RT = (u_tau * (delta / nu)) * std::sqrt(C_f / 2.0);

    // Between Eq. 10 and 11
    Real A1 = 3.7 + 1.5 * beta_c;
    Real A2 = std::min(3.0, 19.0 / std::sqrt(RT)) + 7.0;

    // Eq. 9
    Real C3 = 8.8 * std::pow(RT, -0.57);

    // Eq. 11
    Real F1 = 4.76 * std::pow((1.4 / D), 0.75) * (0.375 * A1 - 1.0);

    // Eq. 13 (hard-coded PI value)
    Real PI = 1.56;

    // Eq. 12 - RHS
    Real term1 = (2.82 * D * D) * std::pow((6.13 * std::pow(D, -0.75) + F1), A1);
    Real term2 = (4.2 * (PI / D) + 1.0) * (omega_bar * omega_bar);
    Real numerator = term1 * term2;

    Real denom1 = std::pow((4.76 * std::pow(omega_bar, 0.75) + F1), A1);
    Real denom2 = std::pow(C3 * omega_bar, A2);
    Real denominator = denom1 + denom2;

    Real RHS = numerator / denominator;

    // Eq. 12 - Phi_pp
    Real Phi_pp = ((tau_w * tau_w) * delta_star / Ux) * RHS;

    return Phi_pp;
}



Real calc_S_qq_Amiet_goody(
    const Real Ux,         // Freestream velocity [m/s]
    const Real omega,      // Angular frequency [rad/s]
    const Real rho,        // Air density [kg/m³]
    const Real tau_w,      // Wall shear stress [Pa]
    const Real delta,      // Boundary layer thickness [m]
    const Real delta_star, // Displacement thickness [m]
    const Real theta,      // Momentum thickness [m] (not used)
    const Real dpdx        // Pressure gradient [Pa/m] (not used)
) {

    
    const Real nu = 1.51e-5;  // Kinematic viscosity [m²/s]

    Real u_tau = std::sqrt(tau_w / rho);
    Real omega_bar = omega * delta / Ux;
    Real Cf = tau_w / (0.5 * rho * Ux * Ux);
    Real RT = (u_tau * delta / nu) * std::sqrt(Cf / 2.0);

    const Real C1 = 0.5;
    const Real C2 = 3.0;
    Real C3 = 1.1 * std::pow(RT, -0.57);

    Real numerator = C2 * omega_bar * omega_bar;
    Real denom = std::pow(std::pow(omega_bar, 0.75) + C1, 3.7)
                 + std::pow(C3 * omega_bar, 7.0);

    Real F = numerator / denom;

    Real Sqq = (tau_w * tau_w * delta * F) / Ux;

    return Sqq;
}