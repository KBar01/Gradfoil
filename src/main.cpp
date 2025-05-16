#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include "data_structs.h"
#include "main_func.h"
#include "panel_funcs.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;

bool runCode(bool restart,Real alphad,const Real Re,const Real Ma,const Real (&xcoords)[Ncoords], Real (&ycoords)[Ncoords],Real (&deltaY)[Ncoords],const Real (&statesInit)[RVdimension],const bool (&turbInit)[Ncoords+Nwake]){

    #if DO_BL_GRADIENT
    Real outputs[10] ; // 10 if doing all gradients CL CD and BL states for both surfaces
    #else
    Real outputs[2] ; // only doing CL and CD gradients  (and fwd output)
    #endif

    // ------------- Doing Adjoint, register relevant input to track gradients --------------------------
    #ifdef USE_CODIPACK

        Tape& tape = Real::getTape();
        tape.setActive();

        for (int i = 0; i < Ncoords; ++i) {
            tape.registerInput(ycoords[i]);
        }
    
    #endif
    //---------------------------- run calculation ----------------------------------------------------

    Real alpha = (alphad/180)*M_PI;
    
    Oper oper(alpha,Re,Ma);
    Geom geom;
    Real flattenedCoords[2 * Ncoords];
    // Flatten the coordinates (column-major order)
    for (int i = 0; i < Ncoords; ++i) {
        flattenedCoords[2*i] =   xcoords[i];     // X[0][i] -> x[i] (x-coordinates)
        flattenedCoords[2*i+1] = ycoords[i];   // X[1][i] -> x[N + i] (y-coordinates)
    }
    Foil foil(flattenedCoords);
    Isol isol;
    Param param;
    Wake wake;
    Vsol vsol;
    
    static Glob glob;
    
    build_gamma_codi(isol,foil,oper);
    init_thermo(oper,param,geom);
    build_wake(foil,geom,oper,isol,wake);
    stagpoint_find(isol,foil,wake);
    identify_surfaces(isol,vsol);
    set_wake_gap(foil,isol,vsol);
    calc_ue_m(foil,wake,isol,vsol);
    rebuild_ue_m(foil,wake,isol,vsol,false);
    

    #ifdef USE_CODIPACK
    
    //init_boundary_layer(oper,foil,param,isol,vsol,glob);
    for (int i=0;i<RVdimension;++i){
        glob.U[i] = statesInit[i];    
    }

    for (int i=0;i<(Ncoords+Nwake);++i){vsol.turb[i] = turbInit[i];}
    
    
    #else
    // Standard run, dealing with restart or not
    if (restart){

        std::ifstream prevfile("restart.json");
        if (!prevfile) {
            std::cerr << "Failed to open restart.json\n";
            return 1;
        }

        // Parse the JSON
        json j;
        prevfile >> j;

        for (int i=0;i<RVdimension;++i){glob.U[i] = j["states"][i];}
        for (int i=0;i<(Ncoords+Nwake);++i){vsol.turb[i] = j["turb"][i];}
    }
    else{
        init_boundary_layer(oper,foil,param,isol,vsol,glob);
    }
    #endif
 
    stagpoint_move(isol,glob,foil,wake,vsol);

    bool converged = solve_coupled(oper,foil,wake,param,vsol,isol,glob);

    Post post;
    calc_force(oper,geom,param,isol,foil,glob,post);

    Real topsurf[4],botsurf[4] ;
    // retusn theta,delta*, tau_max, dp/dx on each surface at 95% chord
    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,vsol,topsurf,botsurf);



    # ifndef USE_CODIPACK
    if (converged){
        json restart;
        restart["states"] = glob.U;
        restart["turb"] = vsol.turb;
        std::ofstream restartFile("restart.json");
        restartFile << restart.dump(4);  // pretty print with 4 spaces indentation
        restartFile.close();

        json out;
        out["CL"]  = post.cl;
        out["CD"]  = post.cd;

        #if DO_BL_GRADIENT
        out["upperMomThickness"] = topsurf[0];
        out["upperDispThickness"] = topsurf[1];
        out["upperTauMax"] = topsurf[2];
        out["upperPressureGrad"] = topsurf[3];

        out["lowerMomThickness"] = botsurf[0];
        out["lowerDispThickness"] = botsurf[1];
        out["lowerTauMax"] = botsurf[2];
        out["lowerPressureGrad"] = botsurf[3];
        #endif
        std::ofstream outFile("out.json");
        outFile << out.dump(4);  // pretty print with 4 spaces indentation
        outFile.close();
    }

    # endif

    // ------------------------ Doing Adjoint: register and store gradients ----------------------------------
    #ifdef USE_CODIPACK

        outputs[0] = post.cl;
        outputs[1] = post.cd;

        #if DO_BL_GRADIENT
        int jacobianHeight = 10;
        #else
        int jacobianHeight = 2;
        #endif

        for (int i=0;i<2;++i){
            tape.registerOutput(outputs[i]);
        }

        #if DO_BL_GRADIENT
        for (int i=0;i<4;++i){
            outputs[2+i] = topsurf[i];
            outputs[2+4+i] = botsurf[i];
        }
        #endif
        
        for (int i=0;i<jacobianHeight;++i){
            tape.registerOutput(outputs[i]);
        }
            
        tape.setPassive();

        for (int i=0;i<jacobianHeight;++i){outputs[i].gradient()[i] = 1.0 ;}
        
        tape.evaluate();

        std::vector<std::string> outputNames = {"CL", "CD",
            "thetaUpper", "deltaStarUpper", "tauMaxUpper", "dpdxUpper",
            "thetaLower", "deltaStarLower", "tauMaxLower", "dpdxLower"};
        
        
        codi::Jacobian<double> jacobian(jacobianHeight,Ncoords);
        
        std::vector<std::vector<double>> allGradients ;
        for (int i = 0; i < Ncoords; ++i) {   
            for (int n=0;n<jacobianHeight;++n){
                jacobian(n,i) = ycoords[i].getGradient()[n];
            }
        }

        // Fill allGradients
        for (int out = 0; out < jacobianHeight; ++out) {
            std::vector<double> grad;
            for (int i = 0; i < Ncoords; ++i) {
                grad.push_back(jacobian(out, i));  // or however you compute the gradient
            }
            allGradients.push_back(grad);
        }


        // Create JSON
        json j;
        for (int i = 0; i < allGradients.size(); ++i) {
            
            j["d " + outputNames[i] + " / d ycoords"] = allGradients[i];
        }

        std::ofstream outFile("ad_gradients.json");
        outFile << j.dump(4);  // pretty-print with 4-space indentation
        outFile.close();
        
        tape.reset();

    #endif

    return converged;
};


int main(){

    // Open JSON file
    std::ifstream infile("input.json");
    if (!infile) {
        std::cerr << "Failed to open input.json\n";
        return 1;
    }

    // Parse the JSON
    json j;
    infile >> j;

    Real xcoords[Ncoords], ycoords[Ncoords];
    Real deltaY[Ncoords] = {0};
    for (int i = 0; i < Ncoords; ++i) {
        xcoords[i] = j["xcoords"][i];       // X[0][i] -> x[i] (x-coordinates)
        ycoords[i] = j["ycoords"][i];   // X[1][i] -> x[N + i] (y-coordinates)
    }
   

    // Read scalar values
    Real targetAlphaDeg = j["alpha_degrees"].get<double>();
    Real Re = j["Re"].get<double>();
    Real Ma = j["Ma"].get<double>();

    int doRestart = j["restart"].get<int>();
    
    Real initStates[RVdimension] = {0};
    bool initTurb[Ncoords+Nwake] = {false} ;

    #ifdef USE_CODIPACK
    
    std::ifstream filein("restart.json");
    if (!filein) {
        std::cerr << "Failed to open states.json\n";
        return 1;
    }

    // Parse the JSON
    json jfile;
    filein >> jfile;

    for (int i=0;i<RVdimension;++i){initStates[i] = jfile["states"][i];}
    for (int i=0;i<(Ncoords+Nwake);++i){initTurb[i] = jfile["turb"][i];}

    
    
    #endif
    
    bool converged = runCode(doRestart,targetAlphaDeg,Re,Ma,xcoords,ycoords,deltaY,initStates,initTurb);
    
    return converged;
};