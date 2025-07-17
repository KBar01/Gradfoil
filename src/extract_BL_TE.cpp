#include <iostream>
#include <cmath>
#include "real_type.h"
#include "get_funcs.h"
#include "data_structs.h"
#include "vector_ops.hpp"
#include "main_func.h"


// Performs cubic interpolation at x using 4 surrounding points (xs, ys)
Real cubic_interp(Real x, const Real* xs, const Real* ys) {
    Real result = 0.0;
    for (int i = 0; i < 4; ++i) {
        Real term = ys[i];
        for (int j = 0; j < 4; ++j) {
            if (i != j)
                term *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += term;
    }
    return result;
}

Real linear_interp(Real x, const Real* xs, const Real* ys) {
    return ys[0] + (ys[1] - ys[0]) * (x - xs[0]) / (xs[1] - xs[0]);
}

Real quadratic_interp(Real x, const Real* xs, const Real* ys) {
    Real result = 0.0;
    for (int i = 0; i < 3; ++i) {
        Real term = ys[i];
        for (int j = 0; j < 3; ++j) {
            if (i != j)
                term *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += term;
    }
    return result;
}

Real adaptive_interp(Real x, const Real* xs, const Real* ys, int n) {
    if (n == 4)
        return cubic_interp(x, xs, ys);
    else if (n == 3)
        return quadratic_interp(x, xs, ys);
    else if (n == 2)
        return linear_interp(x, xs, ys);
    else
        return 0.0; // fallback error
}

int find_interp_position(const Real* xcoords, int start, int end, Real x_target) {

    // For top surface : finds node just before sampling pos
    // For bot surface : finds node just before sampling pos
    
    int found_idx = -1;
    bool is_increasing = xcoords[start] < xcoords[start + 1];

    for (int i = start; i < end - 1; ++i) {
        if (is_increasing) {
            if (xcoords[i] <= x_target && x_target < xcoords[i + 1]) {
                found_idx = i;
                break;
            }
        } else {
            if (xcoords[i] >= x_target && x_target > xcoords[i + 1]) {
                found_idx = i+1;
                break;
            }
        }
    }

    return found_idx;
}

void get_nodes(int topFoundIdx,int botFoundIdx, Real x_target, int* topNodeList,int &topNnodes, int* botNodeList,int &botNnodes, const Vsol&vsol){

    // do top surface :
    int topStart = topFoundIdx - 1 ; // Ideal starting position for cubic interp to have 2 nodes either side of sampling position
    topNnodes = 4 ;
    if (vsol.turb[topFoundIdx] == false){
        topNnodes = 0;
    }
    else{
        if (vsol.turb[topStart] == false){topStart += 1; } // original start is not turbulent, shift start down a node
        topNnodes = (Ncoords-1 - topStart) + 1 ;  // how many available nodes to use for interp
        if (topNnodes>4){topNnodes=4;} // limit to 4 nodes for cubic
    }

    for (int i=0;i<topNnodes;++i){topNodeList[i] = topStart+i ;}


    // do bot surface :
    int botStart = botFoundIdx + 1 ; // Ideal starting position for cubic interp to have 2 nodes either side of sampling position
    botNnodes = 4;
    if (vsol.turb[botFoundIdx] == false){
        botNnodes = 0;
    }
    else{
        if (vsol.turb[botStart] == false){botStart -= 1; } // original start is not turbulent, shift start down a node
        int botNnodes = botStart ;  // how many available nodes to use for interp
        if (botNnodes>4){botNnodes=4;} // limit to 4 nodes for cubic
    }

    for (int i=0;i<botNnodes;++i){botNodeList[i] = botStart-i ;}
}



void interp_BL_states(const int* topIdx,const int* botIdx, const int topNnodes, const int botNnodes, const Real x_target, const Real* xcoords, const Real* states, Real* topInterpStates, Real* botInterpStates){



    // Do top surface first 
    Real txs[4] = {
        xcoords[topIdx[0]],
        xcoords[topIdx[1]],
        xcoords[topIdx[2]],
        xcoords[topIdx[3]]
    };

    for (int q = 0; q < 4; ++q) {
        Real tys[4] = {
            states[colMajorIndex(q,topIdx[0],4)],
            states[colMajorIndex(q,topIdx[1],4)],
            states[colMajorIndex(q,topIdx[2],4)],
            states[colMajorIndex(q,topIdx[3],4)]
        };
        topInterpStates[q] = adaptive_interp(x_target,txs,tys,topNnodes);
    }

    // then bottom surface

    Real bxs[4] = {
        xcoords[botIdx[0]],
        xcoords[botIdx[1]],
        xcoords[botIdx[2]],
        xcoords[botIdx[3]]
    };

    for (int q = 0; q < 4; ++q) {
        Real bys[4] = {
            states[colMajorIndex(q,botIdx[0],4)],
            states[colMajorIndex(q,botIdx[1],4)],
            states[colMajorIndex(q,botIdx[2],4)],
            states[colMajorIndex(q,botIdx[3],4)]
        };
        botInterpStates[q] = adaptive_interp(x_target,bxs,bys,botNnodes);
    }
}


// Computes dp/dx at x_target using cubic interpolation + central finite difference
Real interpolate_dpdx(const Real* xcoords, const Real* Cps, const int* nodeIdx, const int nodeN, Real x_target, const Oper&oper,const Geom&geom,const Real Uinf) {
    
    
    Real h = 1e-6 * geom.chord ;  // Small step size for derivative approximation TODO: verfiy step is correct (convergence)

    // Get x positions slightly left and right of target
    Real x_plus = x_target + h;
    Real x_minus = x_target - h;

    Real xs[4] = {
        xcoords[nodeIdx[0]],
        xcoords[nodeIdx[1]],
        xcoords[nodeIdx[2]],
        xcoords[nodeIdx[3]]
    };

    Real ps[4] = {
        Cps[nodeIdx[0]],
        Cps[nodeIdx[1]],
        Cps[nodeIdx[2]],
        Cps[nodeIdx[3]]
    };

    // Interpolate pressure at x+h and x-h
    Real CpPlus  = adaptive_interp(x_plus,xs,ps,nodeN);
    Real CpMinus = adaptive_interp(x_minus,xs,ps,nodeN);

    // Central finite difference
    Real dpdx = ((CpPlus - CpMinus) / (2.0 * h)) * (0.5 * oper.rho * Uinf*Uinf);

    return dpdx;
}


Real interpolate_cf(const Real* xcoords, const Real* states, const int* nodeIdx, const int nodeN, Real x_target, const Vsol&vsol, const Param&param) {
    
    Real xs[4] = {
        xcoords[nodeIdx[0]],
        xcoords[nodeIdx[1]],
        xcoords[nodeIdx[2]],
        xcoords[nodeIdx[3]]
    };


    Real cf[4];

    int indexes[4] = {nodeIdx[0], nodeIdx[1],nodeIdx[2], nodeIdx[3]} ;
    Real cf_U[4] = {0};
    for (int i=0;i<4;++i){
        cf[i] = get_cf(
            states[colMajorIndex(0,indexes[i],4)],
            states[colMajorIndex(1,indexes[i],4)],
            states[colMajorIndex(2,indexes[i],4)],
            states[colMajorIndex(3,indexes[i],4)],
            vsol.turb[indexes[i]],
            false,
            param,
            cf_U
        );
    }
    // Interpolate 
    Real Cf95  = adaptive_interp(x_target,xs,cf,nodeN);

    return Cf95;
}


void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real*Cps, const Oper&oper, const Vsol&vsol, const Param&param,
    Real (&topBLStates)[7],Real (&botBLStates)[7],const Real Uinf, const Geom&geom,const Real x_target) {

    //Real x_target = 0.95;
    
    /* State order: theta, delta*, tau_max, Ue, dpdx, tau_wall, delta 99% thickness*/
    
    // find index of node before sampling position (top and bottom)
    int foundIndexBot = find_interp_position(xcoords,0,30,x_target);
    int foundIndexTop = find_interp_position(xcoords,Ncoords-30, Ncoords-1,x_target);

    // find the indexes to to the interpolation over 
    int topIdx[4] = {0},botIdx[4] = {0}, topN, botN ;
    get_nodes(foundIndexTop,foundIndexBot,x_target,topIdx,topN,botIdx,botN,vsol);

    interp_BL_states(topIdx,botIdx,topN,botN,x_target,xcoords,states,topBLStates,botBLStates);
    
    Real dpdxBot = interpolate_dpdx(xcoords,Cps,botIdx,botN,x_target,oper,geom,Uinf);
    Real dpdxTop = interpolate_dpdx(xcoords,Cps,topIdx,topN,x_target,oper,geom,Uinf);
    
    // ---------------------- bottom surface dimensionals -----------------------------------------
    Real ignore;
    Real UeCorrected = (get_uk(botBLStates[3],param,ignore)) * Uinf;
    botBLStates[3] = UeCorrected;
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxBot = (botBLStates[2] * botBLStates[2]) * (oper.rho * (botBLStates[3]*botBLStates[3])) ;
    
    botBLStates[2] = tauMaxBot;
    botBLStates[4] = dpdxBot;

    // tau_wall = Cf*(rho * ue^2) / 2
    Real cfBot = interpolate_cf(xcoords,states,botIdx,botN,x_target,vsol,param);
    Real tauWallBot = (cfBot/2) * oper.rho * botBLStates[3] * botBLStates[3] ;
    botBLStates[5] = tauWallBot ;
    
    // now get 99% thickness
    Real deltaBot = 0.0;
    Real frictionVel = std::sqrt(tauWallBot/oper.rho) ;
    if (tauWallBot!=0.0){
        deltaBot = botBLStates[0]*(3.15 + 1.72/((botBLStates[1]/botBLStates[0]) - 1)) + botBLStates[1] ;
    }
    botBLStates[6] = deltaBot;

    // --------------------------------- top surface dimensionals ----------------------------------------
    UeCorrected = (get_uk(topBLStates[3],param,ignore)) * Uinf;
    topBLStates[3] = UeCorrected;
    
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxTop = (topBLStates[2] * topBLStates[2]) * (oper.rho * (topBLStates[3]*topBLStates[3])) ;

    topBLStates[2] = tauMaxTop;
    topBLStates[4] = dpdxTop;
    
    // tau_wall = Cf*(rho * ue^2) / 2
    Real cfTop = interpolate_cf(xcoords,states,topIdx,topN,x_target,vsol,param);
    Real tauWallTop = (cfTop/2) * oper.rho * topBLStates[3] * topBLStates[3] ;
    topBLStates[5] = tauWallTop ;

    // now get 99% thickness
    Real deltaTop = 0.0;
    frictionVel = std::sqrt(tauWallTop/oper.rho) ;
    if (tauWallTop!=0.0){
        deltaTop = topBLStates[0]*(3.15 + 1.72/((topBLStates[1]/topBLStates[0]) - 1)) + topBLStates[1] ;
    }
    
    topBLStates[6] = deltaTop;

}



