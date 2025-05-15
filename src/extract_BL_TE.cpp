#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "data_structs.h"
#include "vector_ops.hpp"
#include "main_func.h"


// Performs cubic interpolation at x using 4 surrounding points (xs, ys)
Real cubic_interp(Real x, const Real xs[4], const Real ys[4]) {
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

int find_interp_position(const Real* xcoords, int start, int end, Real x_target) {

    // Prepare 4-point stencil for x around x_target (shared between both interpolations)
    int found_idx = -1;
    bool is_increasing = xcoords[start] < xcoords[start + 1];

    for (int i = start + 1; i < end - 2; ++i) {
        if (is_increasing) {
            if (xcoords[i] <= x_target && x_target < xcoords[i + 1]) {
                found_idx = i;
                break;
            }
        } else {
            if (xcoords[i] >= x_target && x_target > xcoords[i + 1]) {
                found_idx = i;
                break;
            }
        }
    }

    return found_idx;
}



/*
 * Interpolates the 4 flow quantities at x_target on one surface (top or bottom).
 * Inputs:
 *   - xcoords: array of x-coordinates [size N]
 *   - states: flattened 1D array of 4*N flow states (4 values per node)
 *   - start, end: range of indices defining the surface segment (e.g., 0 to N/2)
 *   - x_target: x-position at which to interpolate (e.g., 0.95 for 95% chord)
 * Outputs:
 *   - interpolated_states: array of 4 values, filled with interpolated results
 */
void interpolate_on_surface(const Real* xcoords, const Real* states, int found_idx, Real x_target, Real* interpolated_states, const Vsol&vsol) {
   
    Real xs[4] = {
        xcoords[found_idx - 1],
        xcoords[found_idx],
        xcoords[found_idx + 1],
        xcoords[found_idx + 2]
    };

    int entries[3] = {0,1,3};
    for (int q = 0; q < 3; ++q) {
        Real ys[4] = {
            states[colMajorIndex(entries[q],found_idx-1,4)],
            states[colMajorIndex(entries[q],found_idx  ,4)],
            states[colMajorIndex(entries[q],found_idx+1,4)],
            states[colMajorIndex(entries[q],found_idx+2,4)]
        };
        interpolated_states[entries[q]] = cubic_interp(x_target, xs, ys);
    }

    // Checking Ctau interp for no transition or laminar node inclusion
    bool bottom = true;
    if (found_idx > 100){bottom = false;}

    bool laminarBL = false; 
    Real ctauX[4], cTau[4] ;

    if (bottom){

        if (vsol.turb[found_idx+1] == false){

            interpolated_states[2] = 0; // hasnt transitione before 95%, no turb BL value
            laminarBL = true ;
        }
        else if (vsol.turb[found_idx+2] == false){
            // the last entry on interp grid would be laminar, step grid back to make sure all turb 
            
            ctauX[0] = xcoords[found_idx -2];
            ctauX[1] = xcoords[found_idx -1];
            ctauX[2] = xcoords[found_idx   ];
            ctauX[3] = xcoords[found_idx +1];

            cTau[0] = states[colMajorIndex(2,found_idx-2,4)];
            cTau[1] = states[colMajorIndex(2,found_idx-1,4)];
            cTau[2] = states[colMajorIndex(2,found_idx  ,4)];
            cTau[3] = states[colMajorIndex(2,found_idx+1,4)];
        }
        else{
            
            // use the normal grid
            ctauX[0] = xcoords[found_idx -1];
            ctauX[1] = xcoords[found_idx   ];
            ctauX[2] = xcoords[found_idx +1];
            ctauX[3] = xcoords[found_idx +2];

            cTau[0] = states[colMajorIndex(2,found_idx-1,4)];
            cTau[1] = states[colMajorIndex(2,found_idx  ,4)];
            cTau[2] = states[colMajorIndex(2,found_idx+1,4)];
            cTau[3] = states[colMajorIndex(2,found_idx+2,4)];
        }
    }
    else {

        if (vsol.turb[found_idx] == false){
            
            interpolated_states[2] = 0; // hasnt transitione before 95%, no turb BL value
            laminarBL = true ;
        }
        else if (vsol.turb[found_idx-1] == false){
            
            // first node in grid was laminar, step forward to make all turb
            ctauX[0] = xcoords[found_idx   ];
            ctauX[1] = xcoords[found_idx +1];
            ctauX[2] = xcoords[found_idx +2];
            ctauX[3] = xcoords[found_idx +3];

            cTau[0] = states[colMajorIndex(2,found_idx  ,4)];
            cTau[1] = states[colMajorIndex(2,found_idx+1,4)];
            cTau[2] = states[colMajorIndex(2,found_idx+2,4)];
            cTau[3] = states[colMajorIndex(2,found_idx+3,4)];

        }
        else{
            // use the normal grid
            ctauX[0] = xcoords[found_idx -1];
            ctauX[1] = xcoords[found_idx   ];
            ctauX[2] = xcoords[found_idx +1];
            ctauX[3] = xcoords[found_idx +2];

            cTau[0] = states[colMajorIndex(2,found_idx-1,4)];
            cTau[1] = states[colMajorIndex(2,found_idx  ,4)];
            cTau[2] = states[colMajorIndex(2,found_idx+1,4)];
            cTau[3] = states[colMajorIndex(2,found_idx+2,4)];
        }
    }

    if (!laminarBL){
        interpolated_states[2] = cubic_interp(x_target, ctauX, cTau);
    }
}


// Computes dp/dx at x_target using cubic interpolation + central finite difference
Real interpolate_dpdx(const Real* xcoords, const Real* Cps, int found_idx, Real x_target, const Oper&oper) {
    
    
    Real h = 1e-6;  // Small step size for derivative approximation TODO: verfiy step is correct (convergence)

    // Get x positions slightly left and right of target
    Real x_plus = x_target + h;
    Real x_minus = x_target - h;

    Real xs[4] = {
        xcoords[found_idx - 1],
        xcoords[found_idx],
        xcoords[found_idx + 1],
        xcoords[found_idx + 2]
    };

    Real ps[4] = {
        Cps[found_idx - 1],
        Cps[found_idx],
        Cps[found_idx + 1],
        Cps[found_idx + 2]
    };

    // Interpolate pressure at x+h and x-h
    Real CpPlus  = cubic_interp(x_plus, xs, ps);
    Real CpMinus = cubic_interp(x_minus, xs, ps);

    // Central finite difference
    Real dpdx = ((CpPlus - CpMinus) / (2.0 * h)) / (0.5 * oper.rho * oper.Vinf*oper.Vinf);

    return dpdx;
}



/*
 * Master function to interpolate 4 flow quantities at 95% chord on both surfaces.
 * Inputs:
 *   - xcoords: array of x-coordinates, size N
 *   - states: flattened array of 4*N states, with 4 values per node
 *   - N: total number of surface nodes (must be even, with bottom and top halves)
 * Outputs:
 *   - interp_bottom: interpolated 4 values on bottom surface at 95% chord
 *   - interp_top: interpolated 4 values on top surface at 95% chord
 */
void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real*Cps, const Oper&oper, const Vsol&vsol,
    Real (&topBLStates)[4],Real (&botBLStates)[4]) {

    Real x_target = 0.95;
    
    //  --------------------------bot surface BL states interp -----------------------------------------
    int foundIndexBot = find_interp_position(xcoords,0,20,x_target);
    interpolate_on_surface(xcoords, states, foundIndexBot, x_target, botBLStates,vsol);
    Real dpdxBot = interpolate_dpdx(xcoords,Cps,foundIndexBot,x_target,oper);
    
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxBot = (botBLStates[2] * botBLStates[2]) * (oper.rho * (botBLStates[3]*botBLStates[3])) ;
    
    botBLStates[2] = tauMaxBot;
    botBLStates[3] = dpdxBot;
    
    // ---------------------------- top surface BL interp ------------------------------------------
    int foundIndexTop = find_interp_position(xcoords,Ncoords-20, Ncoords-1,x_target);
    interpolate_on_surface(xcoords, states,foundIndexTop, x_target, topBLStates,vsol);
    Real dpdxTop = interpolate_dpdx(xcoords,Cps,foundIndexTop,x_target,oper);
    
    // BLstate is C_tau ^ 0.5 , and C_tau = tau_max / (rho * Ue^2)
    Real tauMaxTop = (topBLStates[2] * topBLStates[2]) * (oper.rho * (topBLStates[3]*topBLStates[3])) ;

    topBLStates[2] = tauMaxTop;
    topBLStates[3] = dpdxTop;
}



