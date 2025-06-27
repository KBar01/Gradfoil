#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "residuals.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"






void stagnation_state(const Real*U1,const Real*U2,const Real x1,const Real x2,
    Real (&Ust)[4],Real (&Ust_U)[32],Real (&Ust_x)[8],Real&xst){


    Real dx = x2-x1;
    Real dx_x[2] = {-1, 1};
    Real rx = x2/x1;
    Real rx_x[2] = {-rx/x1,1/x1};
  
    // linear extrapolation weights and stagnation state
    Real w1 =  x2/dx, w1_x[2] = {-w1/dx*dx_x[0], -w1/dx*dx_x[1] + 1/dx};
    Real w2 = -x1/dx, w2_x[2] = {-w2/dx*dx_x[0] -1/dx, -w2/dx*dx_x[1]};
        
    for (int i=0;i<4;++i){Ust[i] = U1[i]*w1 + U2[i]*w2;}
    // quadratic extrapolation of the edge velocity for better slope, ue=K*x
    Real wk1 = rx/dx,       wk1_x[2] = {rx_x[0]/dx - wk1/dx*dx_x[0], rx_x[1]/dx - wk1/dx*dx_x[1]};
    Real wk2 = -1/(rx*dx),  wk2_x[2] = {-wk2*(rx_x[0]/rx + dx_x[0]/dx), -wk2*(rx_x[1]/rx + dx_x[1]/dx)}; 
    Real K = wk1*U1[3] + wk2*U2[3] ;
    Real K_U[8] = {0,0,0,wk1, 0,0,0,wk2};
    Real K_x[2] = {U1[3]*wk1_x[0] + U2[3]*wk2_x[0], U1[3]*wk1_x[1] + U2[3]*wk2_x[1]};
  
    //stagnation coord cannot be zero, but must be small
    xst = 1e-6;
    Ust[3] = K*xst ; // linear dep of ue on x near stagnation
    //Ust_U = np.block([[w1*np.eye(3,4), w2*np.eye(3,4)], [K_U*xst]])
    
    // ---- columns 0‑3 :  w1 * eye(3,4)  ----
    for (int c = 0; c < 4; ++c) {          // 4 columns
        for (int r = 0; r < 4; ++r) {      // 4 rows
            // identity on rows 0‑2 only
            Real val = (r == c && r < 3) ? w1 : 0.0;
            Ust_U[r + c*4] = val;
        }
    }
    // ---- columns 4‑7 :  w2 * eye(3,4)  ----
    for (int c = 0; c < 4; ++c) {
        for (int r = 0; r < 4; ++r) {
            Real val = (r == c && r < 3) ? w2 : 0.0;
            Ust_U[r + (c+4)*4] = val;
        }
    }
    // ---- last row (row 3) : K_U * xst  ----
    for (int c = 0; c < 8; ++c) {
        Ust_U[3 + c*4] = K_U[c] * xst;     // K_U is length‑8 array
    }

    Real temp1[6] = {0}, temp2[6] = {0};
    cnp::outer_product<3,2>(U1,w1_x,temp1);
    cnp::outer_product<3,2>(U2,w2_x,temp2);
    cnp::add_inplace<6>(temp1,temp2);
    
    Real botRow[2];
    cnp::scalar_mul<2>(K_x,xst,botRow);
    cnp::vstack<3,1,2>(temp1,botRow,Ust_x);
    //Ust_x = np.vstack((np.outer(U1[0:3],w1_x) + np.outer(U2[0:3],w2_x), K_x*xst))
 
}


void build_glob_RV(const Foil&foil, const Vsol&vsol,const Isol&isol,Glob&glob, Param&param, Trans&tdata,const Real topNcrit, const Real botNcrit){
    
    constexpr int RVsize = 4*(Ncoords+Nwake);
    constexpr int RXsize = 3*(Ncoords+Nwake);
    
    Real R_st[3*(Ncoords+Nwake)] = {0};


    const Real* xi = isol.distFromStag;
    for (int si = 0; si < 3; ++si) {    // for each surface (upper/lower/wake)
        
        
        if (si != 1){
            param.ncrit = botNcrit;
        }
        else{
            param.ncrit = topNcrit;
        }
        
        const std::vector<int>& Is = vsol.Is[si]; // list of surface node indices from stag point
        const int nSurfPoints = Is.size();

        // Check for edge case of first node hitting stag point exactly
        int i0 = ((si < 2) && (xi[Is[0]] < 1e-8 * xi[Is[nSurfPoints-1]])) ? 1 : 0;

        bool turb = false, wake=false ; 
        
        if (si < 2) {
            
            Real R1[3], R1_U[24], R1_x[6]; // temp storage
            
            /*
            if (3*Is[i0] == 280 || 3*Is[i0]+1 == 280 || 3*Is[i0]+2 == 280 ){
                int nothing_value = 2 + 3+ Is[i0] ;
            }*/

            Real* U1 = &glob.U[colMajorIndex(0,Is[i0],4)];
            Real* U2 = &glob.U[colMajorIndex(0,Is[i0+1],4)];
            
            // Compute stagnation state
            Real Ust[4], Ust_U[32], Ust_x[8], xst; 
            stagnation_state(U1,U2,xi[Is[i0]],xi[Is[i0+1]],Ust,Ust_U,Ust_x,xst);

            Real R1_Ut[24];
            residual_station(Ust,Ust,xst,xst,0,0,false,false,true,param,R1,R1_Ut,R1_x);

            Real R1_Ust[12];
            cnp::add<12>(R1_Ut,R1_Ut+12,R1_Ust);
            
            cnp::matmat_mul<3,4,8>(R1_Ust,Ust_U,R1_U);  // R1_U = R1_Ust @ Ust_U
            cnp::matmat_mul<3,4,2>(R1_Ust,Ust_x,R1_x);  // R1_x = R1_Ust @ Ust_x

            int J[2] = {Is[i0],Is[i0+1]};

            if (i0 == 1) {

                // i0=0 point landed right on stagnation: set value to Ust
                int Ig = 3 * Is[0];
                for (int i = 0; i < 3; ++i)
                    {glob.R[Ig + i] = glob.U[colMajorIndex(i,Is[0],4)] - Ust[i];
                }

                for (int col = 0; col < 4; ++col){
                    for (int row = 0; row < 3; ++row){
                        glob.R_V[colMajorIndex(Ig+row,4*Is[0]+col,RVsize)] += (row == col ? 1.0 : 0.0);
                        glob.R_V[colMajorIndex(Ig+row,4*J[0]+col,RVsize)]  -= Ust_U[colMajorIndex(row,col,4)];
                        glob.R_V[colMajorIndex(Ig+row,4*J[1]+col,RVsize)]  -= Ust_U[colMajorIndex(row,col+4,4)];
                    }
                }
                
                
                // Swap out for building Rs_t directly:
                for (int row=0;row<3;++row){
                    
                    Real firstVal = -Ust_x[colMajorIndex(row,0,4)];
                    if (isol.edgeVelSign[J[0]] == -1){R_st[Ig+row] += firstVal;}
                    else{R_st[Ig+row] -= firstVal;}

                    Real secondVal = -Ust_x[colMajorIndex(row,1,4)];
                    if (isol.edgeVelSign[J[1]] == -1){R_st[Ig+row] += secondVal;}
                    else{R_st[Ig+row] -= secondVal;}
                }
            }
            
            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row){glob.R[Ig+row] = R1[row];}

            // here as well, swapping out the Rx build to build R_st directly
            for (int j = 0; j < 2; ++j){
                
                cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*J[j],R1_U,3,0,4*j,3,4);
                //cnp::equate_block_inplace(glob.R_x,RXsize,Ig,J[j],R1_x,3,0,j,3,1); //TODO seperate out into two things as J not ordered
            }

            for (int row=0;row<3;++row){
                    
                Real firstVal = R1_x[colMajorIndex(row,0,3)];
                if (isol.edgeVelSign[J[0]] == -1){R_st[Ig+row] += firstVal;}
                else{R_st[Ig+row] -= firstVal;}

                Real secondVal = R1_x[colMajorIndex(row,1,3)];
                if (isol.edgeVelSign[J[1]] == -1){R_st[Ig+row] += secondVal;}
                else{R_st[Ig+row] -= secondVal;}
            }
        } 
        else {  // dealing with start of wake
            Real R1[3], R1_U[36]={0};
            int J[3]; 
            wake_sys(vsol,foil,glob,param,R1,R1_U,J);

            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row){glob.R[Ig+row] = R1[row];}

            for (int j = 0; j < 3; ++j){
                cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*J[j],R1_U,3,0,4*j,3,4);
            }

            wake = true;
            turb = true;
        }

        // loop over remaining points in surface
        for (int i = i0 + 1; i < nSurfPoints; ++i) {
            
            int prevI = i - 1;
            int currI = i;
            bool tran = vsol.turb[Is[prevI]] ^ vsol.turb[Is[currI]];
            
            Real Ri[3], Ri_U[3*8], Ri_x[3*2];
            
            Real* Uprev = &glob.U[colMajorIndex(0,Is[prevI],4)];
            Real* Ucurr = &glob.U[colMajorIndex(0,Is[currI],4)];

            if (tran){

                int isForced = tdata.isForced[si];
                Real transPos = tdata.transPos[si];

                // do linear interp to find y value
                Real x1 = foil.x[colMajorIndex(0,Is[prevI],2)],y1 = foil.x[colMajorIndex(1,Is[prevI],2)];
                Real x2 = foil.x[colMajorIndex(0,Is[currI],2)],y2 = foil.x[colMajorIndex(1,Is[currI],2)];

                Real yt = y1 + ((tdata.transPos[si] - x1) / (x2 - x1)) * (y2 - y1);
                
                Real distFromStagTrans =  isol.distFromStag[Is[prevI]] + std::sqrt((tdata.transPos[si]-x1)*(tdata.transPos[si]-x1)  + (yt-y1)*(yt-y1));


                residual_transition(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],0,0,param,isForced,distFromStagTrans,Ri,Ri_U,Ri_x);}
            else{
                Real aux1=0,aux2=0;
                if (wake){aux1=vsol.wgap[Is[prevI]-Ncoords];aux2 = vsol.wgap[Is[currI]-Ncoords] ;}
                residual_station(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],aux1,aux2,wake,turb,false,param,Ri,Ri_U,Ri_x);
            }
            
            // update residuals
            int Ig = 3*Is[i];
            for (int j = 0; j < 3; ++j){glob.R[Ig + j] += Ri[j];}
            
            // update R_U (or R_V in this case)
            cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*Is[i-1],Ri_U,3,0,0,3,4);
            cnp::equate_block_inplace(glob.R_V,RVsize,Ig,4*Is[i  ],Ri_U,3,0,4,3,4);

            // Swap out for building Rs_t directly instead of R_x build first:
            for (int row=0;row<3;++row){
                
                Real firstVal = Ri_x[colMajorIndex(row,0,3)];
                if ( wake || isol.edgeVelSign[Is[i-1]] == 1){R_st[Ig+row] -= firstVal;}
                else{R_st[Ig+row] += firstVal;}

                Real secondVal = Ri_x[colMajorIndex(row,1,3)];
                if (wake || isol.edgeVelSign[Is[i  ]] == 1 ){R_st[Ig+row] -= secondVal;}
                else{R_st[Ig+row] += secondVal;}
            }

            if (tran) {turb = true;}
        }
    }

    // Apply R_x → R_V correction

    cnp::scalar_mul_inplace<RXsize>(R_st,isol.sstag_ue[0]);

    Real* rvColPointer = &glob.R_V[colMajorIndex(0,(4*isol.stagIndex[0] + 3),RVsize)];
    cnp::add_inplace<RXsize>(rvColPointer,R_st);

    Real scale = isol.sstag_ue[1] / isol.sstag_ue[0] ;

    cnp::scalar_mul_inplace<RXsize>(R_st,scale);

    rvColPointer = &glob.R_V[colMajorIndex(0,(4*isol.stagIndex[1] + 3),RVsize)];
    cnp::add_inplace<RXsize>(rvColPointer,R_st);

}


