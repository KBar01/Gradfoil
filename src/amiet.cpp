#include <cmath>
#include <complex>

#include "Faddeeva.hh"
#include "real_type.h"


using std::complex;

using std::sqrt;
using std::abs;
using std::exp;


// Compute sqrt(x + i*y) = u + i*v
void sqrt_complex(Real x, Real y, Real& u, Real& v) {
    Real magnitude = std::sqrt(x * x + y * y);

    // Avoid division by zero or sqrt of negative when magnitude == -x
    Real real_part = std::sqrt((magnitude + x) / 2.0);
    Real imag_part = std::sqrt(std::abs((magnitude - x) / 2.0));

    // Determine sign of imag_part based on input y
    if (y < 0.0) {
        imag_part = -imag_part;
    }

    u = real_part;
    v = imag_part;
}

void compute_Psi_L(
    Real a_r, Real a_i,             // A = a_r + i a_i
    Real s_r, Real s_i,             // sqrt_B_over_B_minus_A = s_r + i s_i
    Real e1_r, Real e1_i,           // erf1 = e1_r + i e1_i
    Real e2_r, Real e2_i,           // erf2 = e2_r + i e2_i
    Real& psi_r, Real& psi_i        // output: real and imag parts of Psi_L
) {
    // === Step 1: Compute exp(2iA) ===
    Real exp2a_mag = std::exp(-2.0 * a_i);
    Real exp2a_r = exp2a_mag * std::cos(2.0 * a_r);
    Real exp2a_i = exp2a_mag * std::sin(2.0 * a_r);

    // === Step 2: Compute 1 / (iA) = (a_i + i a_r) / (a_r^2 + a_i^2) ===
    Real a_mag2 = a_r * a_r + a_i * a_i;
    Real invA_r = a_i / a_mag2;
    Real invA_i = a_r / a_mag2;

    // === Step 3: Compute scalar factor: (-exp(2iA)) / (iA) = i * exp(2iA) / A
    Real tmp_r = exp2a_r * invA_r - exp2a_i * invA_i;
    Real tmp_i = exp2a_r * invA_i + exp2a_i * invA_r;
    Real S_r = -tmp_i;
    Real S_i = tmp_r;

    // === Step 4: Compute exp(-2iA) ===
    Real exp_neg_2a_mag = std::exp(2.0 * a_i);
    Real expm2a_r = exp_neg_2a_mag * std::cos(2.0 * a_r);
    Real expm2a_i = -exp_neg_2a_mag * std::sin(2.0 * a_r);  // -sin(2a_r)

    // === Step 5: Multiply exp(-2iA) * sqrt_B_over_B_minus_A
    Real t1_r = expm2a_r * s_r - expm2a_i * s_i;
    Real t1_i = expm2a_r * s_i + expm2a_i * s_r;

    // === Step 6: Multiply that result by erf1
    Real t2_r = t1_r * e1_r - t1_i * e1_i;
    Real t2_i = t1_r * e1_i + t1_i * e1_r;

    // === Step 7: Final inner expression T = t2 - erf2 + 1
    Real T_r = t2_r - e2_r + 1.0;
    Real T_i = t2_i - e2_i;

    // === Step 8: Final multiply: Psi_L = S * T
    psi_r = S_r * T_r - S_i * T_i;
    psi_i = S_r * T_i + S_i * T_r;
}

void seper_calc_lift(const Real kx,const Real ks,const Real kc,const Real C,const Real MX,const Real beta,const Real mu, Real & psi_LR,Real & psi_LI) {
    
    //Real beta = sqrt(1.0 - MX * MX);
    //Real mu = Mc * kx / (beta * beta);

    Real kx_bar = kx * (C / 2.0);
    Real kc_bar = kc * (C / 2.0);
    Real mu_bar = mu * C / 2.0;

    // Determine if the flow is subcritical or supercritical
    bool is_subcritical = (ks * ks) < (beta * mu) * (beta * mu);

    //complex<Real> kappa;
    Real kappaReal=0, kappaImag=0 ;
    if (is_subcritical) {
        kappaReal = mu * sqrt(1.0 - (ks * ks) / ((beta * mu) * (beta * mu)));
    } 
    else {
        kappaImag = -1.0 * abs(mu) * sqrt((ks * ks) / ((beta * mu) * (beta * mu)) - 1.0);
    }

    //complex<Real> kappa_bar = kappa * C / 2.0;
    Real kappa_barReal = kappaReal*C/2 ;
    Real kappa_barImag = kappaImag*C/2 ;
    
    
    //complex<Real> A = kx_bar + kc_bar;
    Real A = kx_bar + kc_bar ;
    //complex<Real> B = kx_bar + kappa_bar + MX * mu_bar;
    Real BReal = kx_bar + kappa_barReal + MX*mu_bar ;
    Real BImag = kappa_barImag;

    // Compute the square roots and error functions

    // doing B/ (B-A), use idenity for complex division

    Real a = BReal, b = BImag , c = (BReal - A) , d = BImag ;

    Real fractionR = ((a*c) + (b*d)) / ((c*c) + (d*d)) ;
    Real fractionI = ((b*c) + (a*d)) / ((c*c) + (d*d)) ; 


    //complex<Real> sqrt_B_over_B_minus_A = sqrt(B / (B - A));

    Real sqrt_B_over_B_minus_AReal,sqrt_B_over_B_minus_AImag ;
    sqrt_complex(fractionR,fractionI,sqrt_B_over_B_minus_AReal,sqrt_B_over_B_minus_AImag);


    //complex<Real> sqrt_2i_B_minus_A = sqrt(2.0 * I * (B - A));
    fractionR = -2*(BImag) ;
    fractionI = 2*(BReal - A);

    Real erf1InputR,erf1InputI ;

    sqrt_complex(fractionR,fractionI,erf1InputR,erf1InputI);
    
    Real u,v ;
    # ifdef USE_CODIPACK
    const double x_val = erf1InputR.getValue() ;
    const double y_val = erf1InputI.getValue() ;
    const complex<double> z(x_val,y_val) ;
    const complex<double> w = Faddeeva::erf(z);
    
    const std::complex<double> dw_dz = (2.0 / sqrt(M_PI)) * exp(-z * z); // Derivative

    // Compute real and imag parts of Jacobian
    double du_dx = dw_dz.real();        // ∂Re(w)/∂x
    double du_dy = -dw_dz.imag();       // ∂Re(w)/∂y
    double dv_dx = dw_dz.imag();        // ∂Im(w)/∂x
    double dv_dy = dw_dz.real();        // ∂Im(w)/∂y
    
    // Push statement for u = Re(w)
    codi::StatementPushHelper<codi::RealReverse> ph;
    ph.startPushStatement();
    ph.pushArgument(erf1InputR, du_dx);
    ph.pushArgument(erf1InputI, du_dy);
    ph.endPushStatement(u, w.real());


    // Push statement for v = Im(w)

    codi::StatementPushHelper<codi::RealReverse> phIm;
    phIm.startPushStatement();
    phIm.pushArgument(erf1InputR, dv_dx);
    phIm.pushArgument(erf1InputI, dv_dy);
    phIm.endPushStatement(v, w.imag());

    #else

    complex<double> z(erf1InputR,erf1InputI) ;
    complex<double> w = Faddeeva::erf(z);

    u = w.real();
    v = w.imag();

    #endif
    //complex<Real> sqrt_2i_B = sqrt(2.0 * I * B);
    fractionR = -2.0 * (BImag), fractionI = 2.0 * (BReal);
    
    Real erf2InputR, erf2InputI;
    sqrt_complex(fractionR,fractionI,erf2InputR,erf2InputI);
    Real u2,v2 ;

    #ifdef USE_CODIPACK

    double x_val2 = erf2InputR.getValue() ;
    double y_val2 = erf2InputI.getValue() ;
    complex<double> z2(x_val2,y_val2) ;
    complex<double> w2 = Faddeeva::erf(z2);
    
    std::complex<double> dw2_dz2 = (2.0 / sqrt(M_PI)) * exp(-z2 * z2); // Derivative

    // Compute real and imag parts of Jacobian
    du_dx = dw2_dz2.real();        // ∂Re(w)/∂x
    du_dy = -dw2_dz2.imag();       // ∂Re(w)/∂y
    dv_dx = dw2_dz2.imag();        // ∂Im(w)/∂x
    dv_dy = dw2_dz2.real();        // ∂Im(w)/∂y
    
    // Push statement for u = Re(w)
    codi::StatementPushHelper<codi::RealReverse> ph2;
    ph2.startPushStatement();
    ph2.pushArgument(erf2InputR, du_dx);
    ph2.pushArgument(erf2InputI, du_dy);
    ph2.endPushStatement(u2, w2.real());


    // Push statement for v = Im(w)

    codi::StatementPushHelper<codi::RealReverse> ph2Im;
    ph2Im.startPushStatement();
    ph2Im.pushArgument(erf2InputR, dv_dx);
    ph2Im.pushArgument(erf2InputI, dv_dy);
    ph2Im.endPushStatement(v2, w2.imag());

    #else
    complex<double> z2(erf2InputR,erf2InputI) ;
    complex<double> w2 = Faddeeva::erf(z2);

    u2 = w2.real();
    v2 = w2.imag();

    #endif

    
    compute_Psi_L(A,0.0,sqrt_B_over_B_minus_AReal,sqrt_B_over_B_minus_AImag,u,v,u2,v2,psi_LR,psi_LI);
}


// Translated model_lift function
Real model_lift(const Real kx, const Real ks, const Real kc, const Real C, const Real Mc, const Real MX) {
    
    
    const Real beta = std::sqrt(1.0 - MX * MX);
    const Real mu = Mc * kx / (beta * beta);

    // Evaluate lift in both sub/supercritical regimes
    //complex<Real> lift_supercritical = calc_lift(kx,(1e-3 * beta*mu),kc,C,MX,beta,mu);
    Real lift_supercritR,lift_supercritI ;
    //seper_calc_lift(kx,(1e-3 * beta*mu),kc,C,MX,beta,mu,lift_supercritR,lift_supercritI);
    seper_calc_lift(kx,ks,kc,C,MX,beta,mu,lift_supercritR,lift_supercritI);
    //complex<Real> lift_subcritical = calc_lift(kx,(1e3 * beta*mu),kc,C,MX,beta,mu);
    Real lift_subcritR,lift_subcritI ;
    seper_calc_lift(kx,ks,kc,C,MX,beta,mu,lift_subcritR,lift_subcritI);
    
    const Real arg = std::abs(ks / (beta * mu));

    if (arg < 1.0) {
        return std::sqrt((lift_supercritR*lift_supercritR) + (lift_supercritI*lift_supercritI));
    } else {
        if (arg > 10.0) {
            return std::sqrt((lift_subcritR*lift_subcritR) + (lift_subcritI*lift_subcritI));
        } else {
            Real b = std::sqrt((lift_supercritR*lift_supercritR) + (lift_supercritI*lift_supercritI));
            Real m = std::sqrt((lift_subcritR*lift_subcritR) + (lift_subcritI*lift_subcritI)); - b;
            return m * std::log10(arg) + b;
        }
    }
}

// Translated calc_Spp_Freq
Real calc_Spp_Freq(
    const Real c0, const Real rho0, const Real C,const Real MX, const Real omega,
    const Real X, const Real Y, const Real Z,
    const Real S,const Real Phi_qq_input,const int Order
) {
    const Real beta = sqrt(1.0 - (MX*MX));
    const Real s = sqrt(X*X + beta*beta *(Y*Y + Z*Z));
    const Real Ux = MX * c0;
    const Real Mc = 0.8 * MX;
    const Real Uc = Mc * c0;

    const Real kX = omega / Uc;
    const Real kS = (omega / c0) * (Y / s);
    const Real kC = (omega / (c0 / (beta * beta))) * (MX - X/s);

    const Real ls = 1.5 * Uc / omega;

    //Real ls = (1/kX)*(0.62 / (0.62*0.62 + (kS/kX)*(kS/kX)));
    const Real Phi_qq = (1.0 / M_PI) * ls * Phi_qq_input;

    //Real Phi_qq = ls * Phi_qq_input * (Uc/M_PI) ; 
    const Real lift = model_lift(kX, kS, kC, C, Mc, MX);
    const Real b = C / 2.0;

    Real term1 = (omega * b * Z) / (2.0 * M_PI * c0 * s*s);
    
    //Real term1S = (omega/c0)*(C/2)*(Z/(2*M_PI*s*s));
    const Real Spp_p = (term1 * term1) * 2.0*M_PI * S * pow(abs(lift), 2) * Phi_qq;
    //Real Spp_p = (term1 * term1) * (S/2)* (M_PI/Uc) * pow(abs(lift), 2) * Phi_qq *0.3 ;

    
    return Spp_p;
}