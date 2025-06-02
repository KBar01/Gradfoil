#ifndef AMIET_H
#define AMIET_H

#include "real_type.h"

Real calc_S_qq_Amiet_rozenberg(Real Ux, Real omega, Real rho, Real tau_w,
                                 Real delta, Real delta_star, Real theta, Real dpdx);

Real calc_S_qq_Amiet_goody(
    Real Ux,         // Freestream velocity [m/s]
    Real omega,      // Angular frequency [rad/s]
    Real rho,        // Air density [kg/m³]
    Real tau_w,      // Wall shear stress [Pa]
    Real delta,      // Boundary layer thickness [m]
    Real delta_star, // Displacement thickness [m]
    Real theta,      // Momentum thickness [m] (not used)
    Real dpdx        // Pressure gradient [Pa/m] (not used)
);


Real calc_Spp_Freq(
    Real c0, Real rho0, Real C, Real MX, Real omega,
    Real X, Real Y, Real Z,
    Real S, Real Phi_qq_input, int Order
);

#endif
