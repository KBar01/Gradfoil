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
                    Real (&phiqq)[Nsound])
    {


    Real Ue = edgeVel; // will always be non-zero
    
    if (tauWall < 0.0){
        tauWall *= 1.0;
        if (tauWall < 0.0001) {
            tauWall = 0.0001;
        }
    }

    if (tauMax < 0.0){
        tauMax *= 1.0;
        if (tauMax < 0.0001) {
            tauMax = 0.0001;
        }
    }
    
    Real Cf = tauWall/ (0.5*Ue*Ue*rho);
    Real lambda = std::sqrt(2/Cf);
    Real beta_c = std::max((theta/tauWall)*(dpdx),-1.81);
    Real G = 6.1 * std::sqrt(beta_c+1.81) - 1.7;
    
    Real Pi = 0.227;
    if (beta_c > -0.5){
        Pi = 0.8*std::pow(beta_c+0.5, 0.75);
    }
    
    Real H = 1-G/lambda;
    if (H<0.0){
        H = 0.0;
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

// this code follows that by Lee:
/*
Source Characterization of Turbulent
Boundary Layer Trailing Edge Noise Using an
Improved TNO Model
*/

void mean_velocity_profile(const Real (&y)[NblPoints],
                           Real delta,
                           Real u_t,
                           Real nu,          // fluid kinematic viscosity
                           Real Ue,
                           Real (&U)[NblPoints],
                           Real (&dUdy)[NblPoints])
{
    const Real kappa = 0.41;
    const Real B = 5.5;
    for (int i = 0; i < NblPoints; ++i) {
        
        Real y_plus = y[i] * u_t / nu;

        if (y_plus <= 5.0){
            
            //u_plus = y_plus in viscous sub layer 
            U[i]  = u_t * y_plus ;
            // du+/dy = u_t/nu  so  dU/dy = u_t * (u_t/nu) = u_t^2 / nu
            dUdy[i] = (u_t * u_t) / nu;
        }
        else{

        Real W = 1-std::cos(M_PI*y[i] / delta);
        
        Real u_plus = (1.0 / kappa) * std::log(y_plus) + B +
                        0.5*W*((Ue/u_t) - (1/kappa)*std::log((u_t*delta)/nu)-B);
        
        // streamwise velocity
        U[i] = u_plus * u_t;
        // Derivative dU/dy
        // du+/dy:
        Real duplus_dy = (1.0 / (kappa * y[i])) +
                        0.5 * ((Ue / u_t) - (1.0 / kappa) * std::log((u_t * delta) / nu) - B) *
                         (M_PI / delta) * std::sin(M_PI * y[i] / delta);

        dUdy[i] = u_t * duplus_dy;
        }
    }
}

void Turb_shear_stress(
    

    // calculate turbulent shear stress term u2^2 bar
    const Real (&dUdy)[NblPoints],
    const Real (&l_mix)[NblPoints],
    const int isSuction,
    Real (&u22)[NblPoints])
{
   
    for (int i = 0; i < NblPoints; ++i)
    {
       Real nu_t = l_mix[i]*l_mix[i]*std::sqrt(dUdy[i]*dUdy[i]);
       Real kt = std::sqrt((nu_t*nu_t*dUdy[i]*dUdy[i]) / 0.09);
        
       if (isSuction){
        u22[i] = 0.45*kt;
       }
       else{
        u22[i] = 0.3*kt;
       }
    }
}

void Integral_Length_scale(

    // calc vertical integral length scale L_2
    const Real delta,
    const Real (&y)[NblPoints],
    Real (&L2)[NblPoints],
    Real (&l_mix)[NblPoints])
{
    const Real k = 0.38;

    for (int i = 0; i < NblPoints; ++i)
    {
        l_mix[i] = (0.085 * delta * std::tanh( (k*y[i]) / (0.085*delta))) /
                        std::sqrt( std::pow(1+5.5*(y[i]/delta), 6.0) );
        
        L2[i] = l_mix[i] / 0.41 ;
    }
}



void Energy_density_spectrum(
    
    // calculating phi22 for midspan observer in far-field

    const Real k1,
    const Real (&L2)[NblPoints],
    Real (&phi22)[NblPoints])
{

    Real beta1 = 1.0;
    Real beta3 = 0.75;

    for (int i=0;i<NblPoints;++i){

        Real ke = 0.7468 / (2*L2[i]);
        Real term = ((beta1*k1)/ke) * ((beta1*k1)/ke);
        phi22[i] = (4/(9*M_PI)) * ((beta1*beta3)/(ke*ke)) * (term / std::pow(1+term, (7.0/3.0))) ;
    }
}


void calc_WPS_TNO(
    const Real delta,
    Real tauWall,
    const Real edgeVel,
    const Real (&omega)[Nsound],
    const Real rho,
    const Real nu,
    const int isSuction,
    Real (&phiqq)[Nsound])
{
    // Compute shear velocity and min y from target y+
    if (tauWall < 0.0) tauWall *= -1.0;
    Real u_t = std::sqrt(tauWall / rho);

    const Real yplus_target = 0.8;  
    Real y_min = yplus_target * nu / u_t;

    //Build wall-normal grid with cosine stretching
    Real y[NblPoints];
    const Real y_max = delta;
    for (int i = 0; i < NblPoints; ++i)
    {
        Real eta = static_cast<Real>(i) / static_cast<Real>(NblPoints - 1);
        Real y_stretch = 0.5 * (1.0 - std::cos(M_PI * eta));
        y[i] = y_min + (y_max - y_min) * y_stretch;
    }

    Real Uc = 0.65 * edgeVel;

    Real U[NblPoints], dUdy[NblPoints];
    mean_velocity_profile(y, delta, u_t, nu, edgeVel, U, dUdy);

    Real L2[NblPoints], l_mix[NblPoints];
    Integral_Length_scale(delta, y, L2, l_mix);

    Real u22[NblPoints];
    Turb_shear_stress(dUdy, l_mix, isSuction, u22);

    Real phi22[NblPoints];


    for (int w = 0; w < Nsound; ++w)
    {
        Real k1 = omega[w] / Uc;
        Real k = std::abs(k1);

        Energy_density_spectrum(k1, L2, phi22);

        Real integrand[NblPoints];
        for (int i = 0; i < NblPoints; ++i)
        {
  
            Real val = L2[i] * Uc * (dUdy[i] * dUdy[i]) * (u22[i] / (Uc * Uc));
            val *= phi22[i];
            val *= std::exp(-2.0 * y[i] * k);
            
            integrand[i] = val;
        }

        // trapezoidal integration
        Real integral = 0.0;
        for (int i = 1; i < NblPoints; ++i)
        {
            Real dy_local = y[i] - y[i - 1];
            integral += 0.5 * (integrand[i] + integrand[i - 1]) * dy_local;
        }

        // Step 3: prefactor
        Real kfactor = (k1 * k1) / (k * k);
        Real phi_p = 4.0 * rho * rho * kfactor * integral;

        phiqq[w] = phi_p * 2.0;
    }
}