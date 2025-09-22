#include <cmath>
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
                    Real (&phiqq)[Nsound]){
    
    Real a = 3;
    Real b = 2;
    Real c = 0.75;
    Real d = 0.5;
    Real e = 3.7;
    Real f = 1.1;
    Real g = -0.57;
    Real h = 7;
    Real i = 1;
    Real Ue = edgeVel;
    Real u_t = std::sqrt(tauWall/rho);
    Real Rt= (delta/Ue)/(nu/(u_t*u_t));
    Real SS   = Ue / (tauWall*tauWall*delta);
    Real FS   = delta/Ue ;

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        phiqq[n] = ((a*std::pow(omegaBar,b))/(std::pow(i*std::pow(omegaBar, c) + d, e) + std::pow((f*std::pow(Rt, g)*omegaBar), h))) / SS;
    }

}

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
                    Real (&phiqq)[Nsound]){

    Real Ue = edgeVel;
    Real beta_c = std::max((theta/tauWall)*(dpdx),-0.5);

    
    Real Cf = tauWall/ (0.5*Ue*Ue*rho);
    Real lambda = std::sqrt(2/Cf);
    Real G = 6.1 * std::sqrt(beta_c+1.81) - 1.7;
    Real H = 1-G/lambda;

    Real wakeParam = 0.227;
    Real Pi = 0.227;
    if (beta_c > -0.5){
        Pi = 0.8*std::pow(beta_c+0.5, 0.75);
    }
    Real m = 0.5*std::pow(H/1.31, 0.3);
    Real a = 0.45*(1.75*std::pow(Pi*Pi*beta_c*beta_c, m) + 15);
    Real b = 2;
    Real c = 1.637;
    Real d = 0.27;
    Real e = 2.47;
    Real f = std::pow(1.15, -2.0/7.0);
    Real g = -2/7;
    Real h = 7;
    Real i = 1;
    Real u_t = std::sqrt(tauWall/rho);
    Real Rt = (deltaS*u_t*u_t)/(nu*Ue);

    Real SS   = Ue / (tauWall*tauWall*deltaS);
    Real FS   = deltaS/Ue ;

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        phiqq[n] = ((a*std::pow(omegaBar,b))/(std::pow(i*std::pow(omegaBar, c) + d, e) + std::pow((f*std::pow(Rt, g)*omegaBar), h)))/SS;
    }

}

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
                    Real (&phiqq)[Nsound]){

    Real Ue = edgeVel;
    Ue = 64.6;
    
    Real Delta = delta/deltaS  ;
    
    Real beta_c = std::max((theta/tauWall)*(dpdx),-0.5);
    Real Pi = 0.227;
    if (beta_c > -0.5){
        Pi = 0.8*std::pow(beta_c+0.5, 0.75);
    }
    
    /* roz test 
    Delta = 6.0;
    delta = 0.00142;
    deltaS = 0.00236;
    theta = 0.00157;
    tauWall = 5.43;
    tauMax = 5.43;
    beta_c = 3.51;
    Pi = 1.56;
    */

    Real u_t = std::sqrt(tauWall/rho);
    Real Rt = (deltaS*u_t*u_t)/(nu*Ue);



    Real b = 2; // Done
    Real c = 0.75; // Done
    Real A1 = 3.7 + 1.5*beta_c ; //Done - A1
    Real F1 = 4.76*std::pow((1.4/Delta), 0.75) * (0.375*A1 -1) ; // F1
    
    Real a = (2.82*Delta*Delta*std::pow((6.13*std::pow(Delta,-0.75) + F1), A1))  *  (4.2*(Pi/Delta) + 1); //Done
    Real f = 8.8; //done
    Real g = -0.57; //done
    Real F2 = std::min(3.0,19.0/ std::sqrt(Rt)); // done
    Real i = 4.76; //done

    Real SS   = Ue / (tauMax*tauMax*deltaS);
    Real FS   = deltaS/Ue ;

    Real C3prime = 8.8*std::pow(Rt, -0.57);

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        Real top = (a*std::pow(omegaBar, 2));
        Real bot = std::pow(4.76*std::pow(omegaBar, 0.75) + F1, A1)   +   std::pow((C3prime*omegaBar), F2);
        phiqq[n] = ( top/bot )/SS;
    }

}


/////////////////////////////////// All TNO Funcs /////////////////////////////////////

Real dcpdxc_from_dpdx(Real dpdx, Real rho, Real Uinf, Real chord)
{
    // dpdx = ∂p/∂x in Pa/m
    // returns ∂Cp/∂(x/c)
    return (2.0 * chord / (rho * Uinf * Uinf)) * dpdx;
}

void mean_velocity_profile(const Real (&y)[NblPoints],
                           Real delta,
                           Real u_t,
                           Real nu,          // fluid kinematic viscosity
                           Real chord,       // inputs.chord
                           Real Uinf,
                           Real rho,
                           Real dpdx,      // dcpdxc
                           Real tau_w,       // tau_w
                           Real delta_s,     // delta_s
                           Real (&U)[NblPoints],
                           Real (&dUdy)[NblPoints])
{
    const Real k = 0.38;
    const Real B = 5.0;

    Real dcpdxc = dcpdxc_from_dpdx(dpdx,rho,Uinf,chord);

    // dcpdx from nondimensional pressure gradient
    Real dcpdx = dcpdxc * (1.0 / chord);

    // beta parameter
    Real beta = (delta_s / tau_w) * dcpdx;

    // wake parameter
    Real Pi_w = 0.8 * std::pow(beta + 0.5, 0.75);

    // arrays for y_plus and u_plus
    Real y_plus[NblPoints];
    Real u_plus[NblPoints];

    for (int i = 0; i < NblPoints; ++i) {
        y_plus[i] = y[i] * u_t / nu;

        if (y_plus[i] < 5.0) {
            // inner layer
            u_plus[i] = y_plus[i];
        } else {
            // outer layer
            Real outer = (1.0 / k) * std::log(y_plus[i]) + B +
                         (2.0 * Pi_w / k) *
                         std::pow(std::sin(M_PI * y[i] / (2.0 * delta)), 2.0);
            u_plus[i] = outer;
        }

        // streamwise velocity
        U[i] = u_plus[i] * u_t;
    }

    // derivative dUdy
    for (int i = 0; i < NblPoints - 1; ++i) {
        dUdy[i] = (U[i + 1] - U[i]) / (y[i + 1] - y[i]);
    }
    // last point with backward difference
    dUdy[NblPoints - 1] =
        (U[NblPoints - 1] - U[NblPoints - 2]) / (y[NblPoints - 1] - y[NblPoints - 2]);
}

void velocity_fluctuations(const Real* U,
                           Real Uref,
                           Real* u_x,
                           Real* u_y,
                           int N)
{
    const Real gamma = 64.0;
    const Real a = 0.2909;
    const Real b = -0.2598;

    for (int i = 0; i < N; ++i)
    {
        Real Ui = U[i];
        Real Q  = 1.0 - std::exp(-gamma * (1.0 - Ui / Uref));
        u_x[i]  = Ui * ((a + b * Ui / Uref) * Q);
        u_y[i]  = 0.5 * u_x[i];
    }
}

void Integral_Length_scale(
    Real delta,
    const Real (&y)[NblPoints],
    Real (&gamma_y_vv)[NblPoints],
    Real (&gamma_x_uu)[NblPoints])
{
    const Real k = 0.38;

    for (int i = 0; i < NblPoints; ++i)
    {
        Real l_mix = 0.085 * delta * std::tanh((k / 0.085) * y[i] / delta);
        gamma_y_vv[i] = l_mix / k;
        gamma_x_uu[i] = 2.0 * gamma_y_vv[i]; // isotropic turbulence assumption
    }
}


// Compute velocity spectrum at each BL location and frequency
void velocity_spectrum(
    const Real (&gamma_x_uu)[NblPoints],  // integral length scale at each y
    const Real (&omega)[Nsound],          // angular frequencies
    Real Uinf,                            // free stream velocity (inputs.U)
    Real (&phi_uu)[NblPoints][Nsound],    // output
    Real (&phi_vv)[NblPoints][Nsound],    // output
    Real (&kx)[Nsound])                   // streamwise wave number
{
    const Real beta_x = 1.0;
    const Real beta_z = 0.75;

    // constant convective velocity
    Real Uc = 0.7 * Uinf;

    // precompute kx
    for (int n = 0; n < Nsound; ++n)
        kx[n] = omega[n] / Uc;

    // ratio of gamma functions
    const Real gamma_ratio = std::tgamma(5.0/6.0) / std::tgamma(1.0/3.0);
    const Real coeff_uu = gamma_ratio / (std::sqrt(M_PI) * std::tgamma(1.0/3.0));

    for (int i = 0; i < NblPoints; ++i)
    {
        Real ke = (std::sqrt(M_PI) / gamma_x_uu[i]) * gamma_ratio;

        for (int n = 0; n < Nsound; ++n)
        {
            Real kx_over_ke = (beta_x * kx[n]) / ke;

            // phi_uu
            Real denom = std::pow(1.0 + kx_over_ke * kx_over_ke, 5.0/6.0);
            phi_uu[i][n] = coeff_uu * (beta_x / ke) * (1.0 / denom);

            // phi_vv
            Real denom2 = std::pow(1.0 + kx_over_ke * kx_over_ke, 7.0/3.0);
            phi_vv[i][n] = (4.0 / (9.0 * M_PI))
                         * (beta_x * beta_z / (ke * ke))
                         * (kx_over_ke * kx_over_ke)
                         / denom2;
        }
    }
}

void spanwise_correlation_length(
    const Real (&omega)[Nsound], // angular frequencies
    Real Uc,                     // convective velocity
    Real (&gamma_p_z)[Nsound])   // output
{
    const Real bc = 1.4;
    for (int n = 0; n < Nsound; ++n)
    {
        if (omega[n] > 0.0)
            gamma_p_z[n] = bc * Uc / omega[n];
        else
            gamma_p_z[n] = 0.0; // avoid division by zero
    }
}


void Point_spectrum(
    const Real (&U)[NblPoints],
    const Real (&dUdy)[NblPoints],
    const Real (&u_y)[NblPoints],
    const Real (&gamma_y_vv)[NblPoints],
    const Real (&phi_vv)[NblPoints][Nsound], // phi_vv(y,f)
    const Real (&gamma_p_z)[Nsound],
    Real Uc,                                  // convective velocity
    const Real (&omega)[Nsound],
    const Real (&y)[NblPoints],
    Real rho,
    Real kx[Nsound],                          // kx for each frequency
    Real Uinf,                                // freestream speed (inputs.U)
    Real (&Pi_w)[Nsound])                     // output
{
    // compute local convective velocity profile Uc_y(y)
    Real Uc_y[NblPoints];
    for (int j = 0; j < NblPoints; ++j)
    {
        Uc_y[j] = (Uc / Uinf) * U[j];
    }

    for (int i = 0; i < Nsound; ++i)
    {
        // Build temp(y)
        Real temp[NblPoints];
        for (int j = 0; j < NblPoints; ++j)
        {
            Real expTerm = std::exp(-2.0 * std::abs(kx[i]) * y[j]);
            Real num = gamma_y_vv[j] * Uc_y[j] * dUdy[j] * dUdy[j] * (u_y[j] * u_y[j]);
            Real den = (Uc_y[j] * Uc_y[j]);
            temp[j] = (num / den) * phi_vv[j][i] * expTerm;
        }

        // trapz integration over y
        Real integral = 0.0;
        for (int j = 0; j < NblPoints - 1; ++j)
        {
            Real dy = y[j + 1] - y[j];
            integral += 0.5 * dy * (temp[j] + temp[j + 1]);
        }

        Pi_w[i] = (4.0 * M_PI * rho * rho / gamma_p_z[i]) * integral;
    }
}

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
                    Real (&phiqq)[Nsound]){
    
    // 2. Linspace y array
    Real y[NblPoints];
    Real dy = (delta - 0.0001) / (NblPoints - 1);
    for (int i = 0; i < NblPoints; ++i) {
        y[i] = 0.0001 + i * dy;
    }

    // 3. Arrays for intermediate results
    Real U[NblPoints], dUdy[NblPoints];
    Real u_x[NblPoints], u_y[NblPoints];
    Real gamma_y_vv[NblPoints], gamma_x_uu[NblPoints];
    Real phi_uu[NblPoints][Nsound], phi_vv[NblPoints][Nsound];
    Real kx[Nsound], gamma_p_z[Nsound];

    Real u_t = std::sqrt(tauWall/rho);

    //mean velocity profile 
    mean_velocity_profile(y,delta,u_t,nu,chord,Uinf,rho,dpdx,tauWall,deltaS,U,dUdy);

   
    velocity_fluctuations(U,Uinf,u_x,u_y,NblPoints);

    Integral_Length_scale(delta, y, gamma_y_vv, gamma_x_uu);
    
    velocity_spectrum(gamma_x_uu,omega,Uinf,phi_uu,phi_vv,kx);
    
    Real Uc = 0.7*Uinf ;
    spanwise_correlation_length(omega,Uc,gamma_p_z);

    Real Pi_w[Nsound];
    Point_spectrum(U,dUdy,u_y,gamma_y_vv,phi_vv,gamma_p_z,Uc,omega,y,rho,kx,Uinf,Pi_w);
    
    for (int n=0;n<Nsound;++n){
        phiqq[n] = Pi_w[n]*2.0;
    }

}