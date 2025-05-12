#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"

void clear_RV(Glob&glob, const Isol&isol,const Vsol&vsol, const Foil&foil,const Param&param){


    constexpr int RVsize = 4*(Ncoords+Nwake);
    constexpr int RXsize = 3*(Ncoords+Nwake);
    
    const Real* xi = isol.distFromStag;
    for (int si = 0; si < 3; ++si) {    // for each surface (upper/lower/wake)
        
        const std::vector<int>& Is = vsol.Is[si]; // list of surface node indices from stag point
        const int nSurfPoints = Is.size();

        // Check for edge case of first node hitting stag point exactly
        int i0 = ((si < 2) && (xi[Is[0]] < 1e-8 * xi[Is[nSurfPoints-1]])) ? 1 : 0;

        if (si < 2) {
            
            int J[2] = {Is[i0],Is[i0+1]};

            if (i0 == 1) {

                // i0=0 point landed right on stagnation: set value to Ust
                int Ig = 3 * Is[0];
                
                for (int col = 0; col < 4; ++col){
                    for (int row = 0; row < 3; ++row){
                        glob.R_V[colMajorIndex(Ig+row,4*Is[0]+col,RVsize)] = 0;
                        glob.R_V[colMajorIndex(Ig+row,4*J[0]+col,RVsize)]  = 0;
                        glob.R_V[colMajorIndex(Ig+row,4*J[1]+col,RVsize)]  = 0;
                    }
                }
                // building R_x --> can delete if R_st build works directly
            }
            
            int Ig = 3*Is[i0];
            // here as well, swapping out the Rx build to build R_st directly
            for (int j = 0; j < 2; ++j){
                Real R1_U[12]={0};
                cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*J[j],R1_U,3,0,0,3,4);
            }
        } 
        else {  // dealing with start of wake
            Real R1[3], R1_UDummy[36]={0},R1_U[12]={0};
            int J[3]; 
            wake_sys(vsol,foil,glob,param,R1,R1_UDummy,J);

            int Ig = 3*Is[i0];
            for (int j = 0; j < 3; ++j){
                cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*J[j],R1_U,3,0,0,3,4);
            }
        }

        // loop over remaining points in surface
        for (int i = i0 + 1; i < nSurfPoints; ++i) {
            
            Real Ri_U[12]={0};
            // update residuals
            int Ig = 3*Is[i];
            // update R_U (or R_V in this case)
            cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*Is[i-1],Ri_U,3,0,0,3,4);
            cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*Is[i  ],Ri_U,3,0,0,3,4);
        }
    }

    // Apply R_x → R_V correction
    Real* rvColPointer = &glob.R_V[colMajorIndex(0,(4*isol.stagIndex[0] + 3),RVsize)];
    cnp::scalar_mul_inplace<RXsize>(rvColPointer,0);

    rvColPointer = &glob.R_V[colMajorIndex(0,(4*isol.stagIndex[1] + 3),RVsize)];
    cnp::scalar_mul_inplace<RXsize>(rvColPointer,0);


    // now zero ewverything made in the solve_glob_function:
    constexpr int Nsys = Ncoords+Nwake;
    // all edge velocity indices
    int rowStart = 3*Nsys;
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 3;
        for (int row = 0;row<Nsys;++row){
            glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] = 0;
        }
    }

    //all disp thickness indices
    for (int col=0;col<Nsys;++col){

        int colindex = 4*col + 1;
        for (int row = 0;row<Nsys;++row){
            glob.R_V[colMajorIndex(rowStart+row,colindex,4*Nsys)] = 0;
        }
    }
}