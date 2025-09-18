#ifndef AMIET_H
#define AMIET_H

#include "real_type.h"


void calc_WPS_Goody(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real (&phiqq)[Nsound]
);


void calc_WPS_Kamruzzaman(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real (&phiqq)[Nsound]
);

void calc_WPS_Rozenburg(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real (&phiqq)[Nsound]
);

void calc_WPS_TNO(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real chord,
                    Real (&phiqq)[Nsound]
);

Real calc_Spp_Freq(
    Real c0,  Real rho0, Real C,Real MX, Real omega,
    Real X, Real Y, Real Z,
    Real S, Real Phi_qq_input,int Order
);

#endif
