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



Real TNO_velocity_profile(
    Real y2,     // vertical distance [m]
    Real delta,  // boundary layer thickness [m]
    Real ustar,  // friction velocity [m/s]
    Real Ue,     // edge velocity [m/s]
    Real U0,     // free stream velocity [m/s]
    Real nu,     // kinematic viscosity [m^2/s]
    Real kappa,
    Real B
) {
    Real pi = M_PI;
    Real W = 1.0 - std::cos(pi * y2 / delta);
    Real term1 = (1.0 / kappa) * std::log((ustar * y2) / nu) + B;
    Real log_term = std::log((ustar * delta) / nu);
    Real C = (Ue * U0 / ustar) - (1.0 / kappa) * log_term - B;
    Real term2 = 0.5 * W * C;
    return ustar * (term1 + term2);
}

Real analytic_Derivative(
    Real y2,
    Real delta,
    Real ustar,
    Real Ue,
    Real U0,
    Real nu,
    Real kappa,
    Real B
) {
    Real pi = M_PI;
    Real log_term = std::log((ustar * delta) / nu);
    Real C = (U0 * Ue / ustar) - (1.0 / kappa) * log_term - B;

    Real term1 = ustar / (kappa * y2);
    Real term2 = (ustar * C * pi) / (2.0 * delta) * std::sin(pi * y2 / delta);

    return term1 + term2;
}


void compute_TNO_fields(const Real* y2, const Real delta, const Real ustar, const Real Ue, const Real U0, const Real nu,
                        Real* U1, Real* dU,
                        Real* L2, Real* ke, Real* u22, Real kappa, Real B) {
    

    for (int i = 0; i < NblPoints; ++i) {
        Real y2_ratio = y2[i] / delta;

        U1[i] = TNO_velocity_profile(y2[i], delta, ustar, (Ue/U0), U0,nu,kappa,B);
        dU[i] = analytic_Derivative(y2[i],delta,ustar,(Ue/U0),U0,nu,kappa,B);

        // Mixing length
        Real lm_tmp = 0.085 * delta * std::tanh(kappa * y2[i] / (0.085 * delta));
        lm_tmp /= std::sqrt(1.0 + B * std::pow(y2_ratio, 6));
        Real lm = lm_tmp;

        // Eddy viscosity
        Real nut = lm_tmp * lm_tmp * std::abs(dU[i]);

        // Turbulent kinetic energy
        Real kt = std::sqrt(nut * nut * dU[i] * dU[i] / 0.09);

        // L2 length scale
        L2[i] = lm_tmp / kappa;

        // Eddy dissipation constant
        ke[i] = 0.7468 / (2.0 * L2[i]);

        // u2^2
        u22[i] = (2.0 / 3.0) * kt;
    }
}


Real get_phi22(Real k1, Real k3, Real ke, Real beta1, Real beta3) {
    // Non-dimensionalised wavenumbers
    Real k1e = beta1 * k1 / ke;
    Real k3e = beta3 * k3 / ke;

    // Numerator and denominator
    Real numerator = 4.0 * beta1 * beta3;
    Real denominator = 9.0 * M_PI * ke * ke;
    Real denom_power = std::pow(1.0 + k1e * k1e + k3e * k3e, 7.0 / 3.0);

    // Final result
    return ((numerator / denominator) * ((k1e * k1e + k3e * k3e) / denom_power));
}


void compute_surface_pressure_spectrum(
    const Real* omega,        // [Nfreq]
    const Real* y2,          // [NblPoints]
    const Real* dU,          // [NblPoints]
    const Real* L2,          // [NblPoints]
    const Real* ke,    // [NblPoints]
    const Real* u22,         // [NblPoints]
    Real  (&phi)[Nsound],               // [Nfreq] output
    Real rho,
    Real Uinf,
    Real beta1,
    Real beta3
) {
    
    
    Real Uc = 0.7 * Uinf;

    for (int fi = 0; fi < Nsound; ++fi) {


        Real k1 = omega[fi] / Uc;
        Real k3 = 0.0;
        Real kabs = std::sqrt(k1 * k1 + k3 * k3);

        // Compute integrand and integrate via trapezoidal rule
        Real integrand[NblPoints];
        for (int i = 0; i < NblPoints; ++i) {

            Real phi22 = get_phi22(k1, k3, ke[i], beta1, beta3);
            integrand[i] = L2[i] * u22[i] * dU[i] * dU[i] * phi22 * std::exp(-2.0 * kabs * y2[i]);
        }

        // Manual trapezoidal integration over y2
        Real total = 0.0;
        for (int i = 0; i < NblPoints - 1; ++i) {
            Real dy = y2[i + 1] - y2[i];
            total += 0.5 * dy * (integrand[i] + integrand[i + 1]);
        }

        Real Lambda3 = (Uc / omega[fi]) * 1.4; // b = 1
        phi[fi] = ((4.0 * M_PI * rho * rho) / (Lambda3 * Uc)) * total;
    }
}



void calc_S_qq_Amiet_TNO(const Real (&omega)[Nsound],const Real Ue, const Real Uinf, const Real rho,const Real tau_w,
    const Real delta,const Real delta_star,const Real theta,const Real dpdx, Real (&phi)[Nsound]){

    const Real nu = 1.51e-5;
    const Real ustar = std::sqrt(tau_w/rho);
    //create y2
    Real y2[NblPoints];
    
    Real y_start = 1e-6;
    Real y_end   = delta - 1e-7;
    Real dy = (y_end - y_start) / (NblPoints - 1);

    Real add = 0.0 ;
    for (int i = 0; i < NblPoints; ++i) {
        y2[i] = y_start + add * dy;
        add += 1.0;
    }
        
    Real U1[NblPoints];
    Real dU[NblPoints];
    Real u22[NblPoints];
    Real L2[NblPoints];
    Real ke[NblPoints];

    compute_TNO_fields(y2,delta,ustar,Ue,Uinf,nu,U1,dU,L2,ke,u22,0.4,5.3);
    compute_surface_pressure_spectrum(omega,y2,dU,L2,ke,u22,phi,rho,Uinf,1.0,0.75);

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


