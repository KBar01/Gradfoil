#ifndef AMIET_H
#define AMIET_H

#include "real_type.h"

Real calc_S_qq_Amiet_rozenberg(const Real Ux,const Real omega,const Real rho,const Real tau_w,
    const Real delta,const Real delta_star,const Real theta,const Real dpdx);

Real calc_S_qq_Amiet_goody(
    const Real Ux,         // Freestream velocity [m/s]
    const Real omega,      // Angular frequency [rad/s]
    const Real rho,        // Air density [kg/m³]
    const Real tau_w,      // Wall shear stress [Pa]
    const Real delta,      // Boundary layer thickness [m]
    const Real delta_star, // Displacement thickness [m]
    const Real theta,      // Momentum thickness [m] (not used)
    const Real dpdx        // Pressure gradient [Pa/m] (not used)
);


Real calc_Spp_Freq(
    const Real c0, const Real rho0, const Real C,const Real MX, const Real omega,
    const Real X, const Real Y, const Real Z,
    const Real S,const Real Phi_qq_input,const int Order
);

#endif
