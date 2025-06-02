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

bool runCode(bool restart,bool xfoilStart,bool doGetPoints,Real alphad, Real Re, Real Ma,const Real (&inXcoords)[Nin], Real (&inYcoords)[Nin], const Real (&statesInit)[RVdimension],const bool (&turbInit)[Ncoords+Nwake],const Real sampleTE){

    auto start = std::chrono::high_resolution_clock::now();
    #if DO_BL_GRADIENT
    Real outputs[16] ; // 12 if doing all gradients CL CD and BL states for both surfaces
    #else
    Real outputs[2] ; // only doing CL and CD gradients  (and fwd output)
    #endif

    // ------------- Doing Adjoint, register relevant input to track gradients --------------------------
    #ifdef USE_CODIPACK

        Tape& tape = Real::getTape();
        tape.setActive();

        for (int i = 0; i < Nin; ++i) {
            tape.registerInput(inYcoords[i]);
        }
        tape.registerInput(alphad);
        tape.registerInput(Re);
        tape.registerInput(Ma);
    #endif
    //---------------------------- run calculation ----------------------------------------------------

    Real alpha = (alphad/180)*M_PI;
    
    const Real rhoInf = 1.225;
    const Real dynViscInf = 1.789e-5 ;
    
    Real minX = 0.5, maxX = 0.5;
    for (int i=0;i<Nin;++i){
        Real newMin = std::min(minX,inXcoords[i]);
        Real newMax = std::max(maxX,inXcoords[i]);
        minX = newMin ;
        maxX = newMax ;
    }
    Real chordScale = maxX - minX ;
    
    
    Oper oper(alpha,Re,Ma);
    Real Uinf = (Re*dynViscInf)/(oper.rho*chordScale) ; // for scaling the BL outputs later

    Geom geom;
    geom.chord = chordScale;
    
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    Real flattenedCoords[2 * Ncoords]={0};
    make_panels(inCoords,flattenedCoords); // does spline to redist nodes over aerofoil for fixed number of 200 nodes

     
    
    if (doGetPoints){
        
        #ifndef USE_CODIPACK
        json points;
        points["points"]  = flattenedCoords;
        std::ofstream pointsFile("innerfoilNodes.json");
        pointsFile << points.dump(4);  // pretty print with 4 spaces indentation
        pointsFile.close();
        #endif
        return true;
    }
    else
    {
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
    
    else if (xfoilStart){
        init_boundary_layer_from_xfoil(oper,foil,param,isol,vsol,glob);
    }
    
    else {
        init_boundary_layer(oper,foil,param,isol,vsol,glob);
    }
    #endif
 
    stagpoint_move(isol,glob,foil,wake,vsol);

    bool converged = solve_coupled(oper,foil,wake,param,vsol,isol,glob);

    Post post;
    calc_force(oper,geom,param,isol,foil,glob,post);

    Real topsurf[7],botsurf[7];
    // retusn theta,delta*, tau_max, dp/dx, ue on each surface at 95% chord

    Real xcoords[Ncoords]={0};
    Real ycoords[Ncoords]={0};

    for (int i=0;i<Ncoords;++i){
        xcoords[i] = flattenedCoords[colMajorIndex(0,i,2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1,i,2)];
    }

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,vsol,param,topsurf,botsurf,Uinf,geom,(sampleTE*geom.chord));

    auto end = std::chrono::high_resolution_clock::now();

    // Duration in milliseconds (or other units)
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Elapsed time: " << duration.count() << " ms\n";


    std::vector<std::string> outputNames = {"CL", "CD",
        "thetaUpper", "deltaStarUpper", "tauMaxUpper","edgeVelocityUpper", "dpdxUpper", "tauWallUpper", "delta99Upper",
        "thetaLower", "deltaStarLower", "tauMaxLower","edgeVelocityLower", "dpdxLower", "tauWallLower", "delta99Lower"
    };

    # ifndef USE_CODIPACK
    if (converged){
        json restart;
        restart["states"] = glob.U;
        restart["turb"] = vsol.turb;
        std::ofstream restartFile("restart.json");
        restartFile << restart.dump(4);  // pretty print with 4 spaces indentation
        restartFile.close();

        json out;
        out[outputNames[0]]  = post.cl;
        out[outputNames[1]]  = post.cd;

        #if DO_BL_GRADIENT
        out[outputNames[2]] = topsurf[0];
        out[outputNames[3]] = topsurf[1];
        out[outputNames[4]] = topsurf[2];
        out[outputNames[5]] = topsurf[3];
        out[outputNames[6]] = topsurf[4];
        out[outputNames[7]] = topsurf[5];
        out[outputNames[8]] = topsurf[6];

        out[outputNames[9]] = botsurf[0];
        out[outputNames[10]] = botsurf[1];
        out[outputNames[11]] = botsurf[2];
        out[outputNames[12]] = botsurf[3];
        out[outputNames[13]] = botsurf[4];
        out[outputNames[14]] = botsurf[5];
        out[outputNames[15]] = botsurf[6];
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
        constexpr int jacobianHeight = 16;
        #else
        constexpr int jacobianHeight = 2;
        #endif

        #if DO_BL_GRADIENT
        for (int i=0;i<7;++i){
            outputs[2+i] = topsurf[i];
            outputs[2+7+i] = botsurf[i];
        }
        #endif
        
        for (int i=0;i<jacobianHeight;++i){
            tape.registerOutput(outputs[i]);
        }
            
        tape.setPassive();

        for (int i=0;i<jacobianHeight;++i){outputs[i].gradient()[i] = 1.0 ;}
        
        tape.evaluate();

        codi::Jacobian<double> jacobian(jacobianHeight,Nin);
        codi::Jacobian<double> jacobianAlpha(jacobianHeight,1);
        codi::Jacobian<double> jacobianRe(jacobianHeight,1);
        codi::Jacobian<double> jacobianMa(jacobianHeight,1);
        
        std::vector<std::vector<double>> allGradients ;
        for (int i = 0; i < Nin; ++i) {   
            for (int n=0;n<jacobianHeight;++n){
                jacobian(n,i) = inYcoords[i].getGradient()[n];
            }
        }

        for (int n=0;n<jacobianHeight;++n){
            jacobianAlpha(n,0) = alphad.getGradient()[n];
            jacobianRe(n,0)    = Re.getGradient()[n];
            jacobianMa(n,0)    = Ma.getGradient()[n];
        }


        // Fill allGradients
        for (int out = 0; out < jacobianHeight; ++out) {
            std::vector<double> grad;
            for (int i = 0; i < Nin; ++i) {
                grad.push_back(jacobian(out, i));  // or however you compute the gradient
            }
            allGradients.push_back(grad);
        }
        // Create JSON
        json j;
        for (int i = 0; i < allGradients.size(); ++i) {
            j["d " + outputNames[i] + " / d ycoords"] = allGradients[i];
        }
        
        std::vector<double> allGradientsAlf ;
        std::vector<double> allGradientsRe ;
        std::vector<double> allGradientsMa ;
        for (int out = 0; out<jacobianHeight; ++out) {
            allGradientsAlf.push_back(jacobianAlpha(out, 0));
            allGradientsRe.push_back(jacobianRe(out, 0));
            allGradientsMa.push_back(jacobianMa(out, 0));
        }
        
        for (int i = 0; i < allGradientsAlf.size(); ++i) {
            j["d " + outputNames[i] + " / d alpha"] = allGradientsAlf[i];
            j["d " + outputNames[i] + " / d alpha"] = allGradientsRe[i];
            j["d " + outputNames[i] + " / d alpha"] = allGradientsMa[i];
        }
        
        std::ofstream outFile("ad_gradients.json");
        outFile << j.dump(4);  // pretty-print with 4-space indentation
        outFile.close();
        
        tape.reset();

    #endif

    return converged;
    }
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

    Real inXcoords[Nin]={0}, inYcoords[Nin]={0};

    for (int i = 0; i < Nin; ++i) {
        inXcoords[i] = j["xcoords"][i];       // X[0][i] -> x[i] (x-coordinates)
        inYcoords[i] = j["ycoords"][i];   // X[1][i] -> x[N + i] (y-coordinates)
    }
   

    // Read scalar values
    Real targetAlphaDeg = j["alpha_degrees"].get<double>();
    Real Re = j["Re"].get<double>();
    Real Ma = j["Ma"].get<double>();

    int doRestart = j["restart"].get<int>();
    int doXfoilStart = j["xfoilstart"].get<int>();
    int doGetPoints = j["xfoilgetpoints"].get<int>();
    Real sampleTE = j["sampleTE"].get<double>();
    
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
    
    bool converged = runCode(doRestart,doXfoilStart,doGetPoints,targetAlphaDeg,Re,Ma,inXcoords,inYcoords,initStates,initTurb,sampleTE);
    
    std::cout << "converged: " << converged << std::endl ;
    return converged;
};