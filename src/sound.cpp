#include <iostream>
#include <cmath>
#include "data_structs.h"
#include "amiet.h"
#include "real_type.h"
#include <fstream>
#include <sstream>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;


void calc_Spp_Rozenburg(const Real theta,const Real deltaStar,const Real delta,const Real tau,const Real edgeVel,const Real dpdx, const Real (&omega)[Nsound],Real (&Spp)[Nsound],const Oper&oper,const Geom&geom,const Real Uinf,const Real X,const Real Y,const Real Z,const  Real S, Real (&phiqq)[Nsound]){

    // start exp : 2 (100Hz)
    // final exp : 4.30103 (20,000 Hz)


    // loop over this with different omega vals
    
    
    for (int i=0;i<Nsound;++i){
        phiqq[i] = calc_S_qq_Amiet_rozenberg(edgeVel,omega[i],oper.rho,tau,delta,deltaStar,theta,dpdx);
    }

    for (int i=0;i<Nsound;++i){
        Spp[i]  = calc_Spp_Freq(340, oper.rho, geom.chord, (Uinf/340), omega[i], X, Y, Z, S, phiqq[i], 0);
    }
}

void calc_Spp_Goody(const Real theta,const Real deltaStar,const Real delta,const Real tau,const Real edgeVel,const Real dpdx, const Real (&omega)[Nsound],Real (&Spp)[Nsound],const Oper&oper,const Geom&geom,const Real Uinf, const Real X,const Real Y,const Real Z, const Real S){

    // start exp : 2 (100Hz)
    // final exp : 4.30103 (20,000 Hz)


    // loop over this with different omega vals
    Real phiqq[Nsound] ;
    
    for (int i=0;i<Nsound;++i){
        phiqq[i] = calc_S_qq_Amiet_goody(Uinf,omega[i],oper.rho,tau,delta,deltaStar,theta,dpdx);
    }

    for (int i=0;i<Nsound;++i){
        Spp[i]  = calc_Spp_Freq(340, oper.rho, geom.chord, (Uinf/340), omega[i], X, Y, Z, S, phiqq[i], 0);
    }
}


Real calc_OASPL(const Real* botStates, const Real* topStates,const Oper&oper,const Geom&geom, const Real Uinf, const Real X,const Real Y,const Real Z, const Real S,const int doCps){

    const Real startExp = 2.0; // start exp : 2 (100Hz)
    const Real endExp = 4.30103; // final exp : 4.30103 (20,000 Hz)

    Real omega[Nsound];
    Real Freq[Nsound];

    for (int i = 0; i < Nsound; ++i) {
        Real exponent = startExp + (endExp - startExp) * i / (Nsound - 1);
        Freq[i]  = std::pow(10.0, exponent);
        omega[i] =  Freq[i] * 2.0 * M_PI;
    }

    Real SppUpper[Nsound]={0}, SppLower[Nsound]={0};
    Real phiqqUpper[Nsound]={0},phiqqLower[Nsound]={0};
    
    Real theta = topStates[0];
    Real deltaS = topStates[1];
    Real tauMax = topStates[2];
    Real edgeVel = topStates[3];
    Real dpdx = topStates[4];
    Real tauWall = topStates[5];
    Real delta = topStates[6];    
    
    if (tauMax > 0.0){ 
    calc_Spp_Rozenburg(theta,deltaS,delta,tauMax,edgeVel,dpdx,omega,SppUpper,oper,geom,Uinf,X,Y,Z,S,phiqqUpper);
    }

    theta = botStates[0];
    deltaS = botStates[1];
    tauMax = botStates[2];
    edgeVel = botStates[3];
    dpdx = botStates[4];
    tauWall = botStates[5];
    delta = botStates[6];  

    if (tauMax > 0.0){ 
    calc_Spp_Rozenburg(theta,deltaS,delta,tauMax,edgeVel,dpdx,omega,SppLower,oper,geom,Uinf,X,Y,Z,S,phiqqLower);
    }
    
    Real SppTotal[Nsound];
    for (int i=0;i<Nsound;++i){
        
        Real SppSum = SppLower[i] + SppUpper[i];
        SppTotal[i] = SppSum*2*M_PI;
    }

    // Trapezoidal integration over frequency
    Real integral = 0.0;
    for (int i = 0; i < Nsound-1; ++i) {
        Real df = Freq[i + 1] - Freq[i];
        integral += 0.5 * (SppTotal[i] + SppTotal[i + 1]) * df;
    }

    // Convert to OASPL (in dB)
    Real pref2 = (20e-6)*(20e-6);
    Real OASPL = 10.0 * std::log10(integral / pref2);
    #ifndef USE_CODIPACK
    if (doCps){
        
        json amiet;
        amiet["omega"] = omega;
        amiet["phiqqupper"] = phiqqUpper;
        amiet["phiqqlower"] = phiqqLower;
        amiet["sppupper"]   = SppUpper;
        amiet["spplower"]   = SppLower;
        std::ofstream amietFile("amiet.json");
        amietFile << amiet.dump(4);  // pretty print with 4 spaces indentation
        amietFile.close();
    }

    #endif
    return OASPL;
}