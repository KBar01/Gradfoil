#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include "data_structs.h"
#include "main_func.h"
#include "get_funcs.h"
#include "panel_funcs.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;

bool runCode(
    bool restart,
    bool xfoilStart,
    bool doGetPoints,
    Real alphad,
    Real Re, 
    Real Ma,
    const Real (&inXcoords)[Nin], 
    Real (&inYcoords)[Nin],
    const Real (&statesInit)[RVdimension],
    const bool (&turbInit)[Ncoords+Nwake],
    const Real sampleTE,
    const Real X,
    const Real Y,
    const Real Z,
    const Real S,
    Real Uinf,
    const int useCustUinf,
    const int doCps,
    const Real nCrit,
    const Real Ufac, 
    const Real TEfac,
    const Real &topTransPos,
    const Real &botTransPos,
    const bool force,
    const bool Ncrit,
    const Real topNcrit,
    const Real botNcrit){

    auto start = std::chrono::high_resolution_clock::now();
    #if DO_BL_GRADIENT
    Real outputs[16] ; // 12 if doing all gradients CL CD and BL states for both surfaces
    #elif DO_SOUND
    Real outputs;
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
    
    Real minX = 0.5, maxX = 0.01;
    for (int i=0;i<Nin;++i){
        Real newMin = std::min(minX,inXcoords[i]);
        Real newMax = std::max(maxX,inXcoords[i]);
        minX = newMin ;
        maxX = newMax ;
    }
    Real chordScale = maxX - minX ;
    
    
    Oper oper(alpha,Re,Ma);
    
    
    if (!useCustUinf){
        Uinf = (Re*dynViscInf)/(oper.rho*chordScale) ; // for scaling the BL outputs later
    }
    Geom geom;
    geom.chord = chordScale;

    Real flattenedCoords[2 * Ncoords]={0};
    
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    make_panels(inCoords,flattenedCoords,Ufac,TEfac); // does spline to redist nodes over aerofoil for fixed number of 200 nodes
    
    // finding node positions to force transition ------------------------
    Real xTransTop = topTransPos * geom.chord ;
    Real xTransBot = botTransPos * geom.chord ;
    
    int idx_closest_bot = 0;
    int idx_closest_top = Ncoords - 1;  // top TE node
    
    if (force) {
        
        // Bottom surface: from bottom TE forward
        for (int i = 1; i < Ncoords; ++i) {
            Real x = flattenedCoords[colMajorIndex(0, i, 2)];
            Real dist = x - xTransBot;
            if (dist < 0.0) {
                idx_closest_bot = i;
                break;
            }
        }

        // Top surface: from top TE backward
        for (int i = Ncoords - 2; i >= 0; --i) {
            Real x = flattenedCoords[colMajorIndex(0, i, 2)];
            Real dist = x - xTransTop;
            if (dist < 0.0) {
                idx_closest_top = i;
                break;
            }
        }
    }

    Trans tdata;

    tdata.transNode[0] = idx_closest_bot;
    tdata.transNode[1] = idx_closest_top;
    tdata.transPos[0] = xTransBot ;
    tdata.transPos[1] = xTransTop ;

    // -------------------------------------------------------------------------

    //if (doGetPoints){
        
        #ifndef USE_CODIPACK
        json points;
        points["points"]  = flattenedCoords;
        std::ofstream pointsFile("innerfoilNodes.json");
        pointsFile << points.dump(4);  // pretty print with 4 spaces indentation
        pointsFile.close();
        #endif
        //return true;
    //}
    //else
    //{
    Foil foil(flattenedCoords);
    Isol isol;
    Param param;
    param.ncrit = nCrit;
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
        
        if (force){
            tdata.isForced[0]  = 1;
            tdata.isForced[1]  = 1;
        }
    }
    
    //else if (xfoilStart){
        //init_boundary_layer_from_xfoil(oper,foil,param,isol,vsol,glob);
    //}
    
    else {
        init_boundary_layer(oper,foil,param,isol,vsol,glob,tdata,force,topNcrit,botNcrit);
    }
    #endif

    stagpoint_move(isol,glob,foil,wake,vsol);

    bool converged = solve_coupled(oper,foil,wake,param,vsol,isol,glob,tdata,force,topNcrit,botNcrit);

    Post post;
    calc_force(oper,geom,param,isol,foil,glob,post);

    #ifndef USE_CODIPACK
    Real Cf_Uinf[Ncoords];
    if (doCps){
        Real cf_U[4]={0};
        for (int i=0;i<Ncoords;++i){
            Cf_Uinf[i] = get_cf(glob.U[colMajorIndex(0,i,4)],
                                glob.U[colMajorIndex(1,i,4)],
                                glob.U[colMajorIndex(2,i,4)],
                                glob.U[colMajorIndex(3,i,4)],
                                vsol.turb[i],
                                false,
                                param,
                                cf_U
            );
            
            Cf_Uinf[i] *= ((glob.U[colMajorIndex(3,i,4)] * glob.U[colMajorIndex(3,i,4)]));
        }
    }
    #endif


    Real topsurf[7],botsurf[7];
    // retusn theta,delta*, tau_max, dp/dx, ue on each surface at 95% chord

    Real xcoords[Ncoords]={0};
    Real ycoords[Ncoords]={0};

    for (int i=0;i<Ncoords;++i){
        xcoords[i] = flattenedCoords[colMajorIndex(0,i,2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1,i,2)];
    }

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,vsol,param,topsurf,botsurf,Uinf,geom,(sampleTE*geom.chord));


    // if codipack, only use sound code if sound flag on. if not codipack run sound regardless
    #ifdef USE_CODIPACK
    #if DO_SOUND
    Real OASPL = calc_OASPL(botsurf,topsurf,oper,geom,Uinf,X,Y,Z,S,doCps);
    #endif
    #else
    Real OASPL = calc_OASPL(botsurf,topsurf,oper,geom,Uinf,X,Y,Z,S,doCps);
    #endif

    auto end = std::chrono::high_resolution_clock::now();

    // Duration in milliseconds (or other units)
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Elapsed time: " << duration.count() << " ms\n";

    #if DO_SOUND
    //std::vector<std::string> outputNames = {"CL", "CD", "OASPL"};
    std::vector<std::string> outputNames = {"OASPL"};
    #else
    std::vector<std::string> outputNames = {"CL", "CD",
        "thetaUpper", "deltaStarUpper", "tauMaxUpper","edgeVelocityUpper", "dpdxUpper", "tauWallUpper", "delta99Upper",
        "thetaLower", "deltaStarLower", "tauMaxLower","edgeVelocityLower", "dpdxLower", "tauWallLower", "delta99Lower"
    };
    #endif
    

    # ifndef USE_CODIPACK
    if (converged){
        json restart;
        restart["states"] = glob.U;
        restart["turb"] = vsol.turb;
        std::ofstream restartFile("restart.json");
        restartFile << restart.dump(4);  // pretty print with 4 spaces indentation
        restartFile.close();

        json out;
        out["aerofoilChord"] = chordScale;
        out["freestreamVelocity"] = Uinf;
        out["CL"]  = post.cl;
        out["CD"]  = post.cd;

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
        #elif DO_SOUND
        out["OASPL"] = OASPL;
        #endif


        if (doCps){

            //  calc transition point
            Real botTransX = geom.chord;
            for(int i=0;i<isol.stagIndex[0];++i){
                int isTurb = vsol.turb[isol.stagIndex[0] - i];
                if (isTurb){
                    botTransX = foil.x[colMajorIndex(0,isol.stagIndex[0]-i,2)];
                    break;
                } 
            }
            Real topTransX = geom.chord;
            for(int i=0;i<200-isol.stagIndex[1];++i){
                int isTurb = vsol.turb[isol.stagIndex[1] + i];
                if (isTurb){
                    topTransX = foil.x[colMajorIndex(0,isol.stagIndex[1]+i,2)];
                    break;
                } 
            }
            out["innerFoil"] = foil.x;
            out["Cp"] = post.cp;
            out["stagnation"] = isol.stagIndex;
            out["topTransX"]  = topTransX;
            out["botTransX"]  = botTransX;
            
            out["Cf"] = Cf_Uinf;
            std::vector<std::string> BLoutputNames = {"CL", "CD",
                "thetaUpper", "deltaStarUpper", "tauMaxUpper","edgeVelocityUpper", "dpdxUpper", "tauWallUpper", "delta99Upper",
                "thetaLower", "deltaStarLower", "tauMaxLower","edgeVelocityLower", "dpdxLower", "tauWallLower", "delta99Lower"
            };
            out[BLoutputNames[2]] = topsurf[0];
            out[BLoutputNames[3]] = topsurf[1];
            out[BLoutputNames[4]] = topsurf[2];
            out[BLoutputNames[5]] = topsurf[3];
            out[BLoutputNames[6]] = topsurf[4];
            out[BLoutputNames[7]] = topsurf[5];
            out[BLoutputNames[8]] = topsurf[6];

            out[BLoutputNames[9]] = botsurf[0];
            out[BLoutputNames[10]] = botsurf[1];
            out[BLoutputNames[11]] = botsurf[2];
            out[BLoutputNames[12]] = botsurf[3];
            out[BLoutputNames[13]] = botsurf[4];
            out[BLoutputNames[14]] = botsurf[5];
            out[BLoutputNames[15]] = botsurf[6];
        }

        std::ofstream outFile("out.json");
        outFile << out.dump(4);  // pretty print with 4 spaces indentation
        outFile.close();

    }

    # endif

    // ------------------------ Doing Adjoint: register and store gradients ----------------------------------
    #ifdef USE_CODIPACK


        #if DO_BL_GRADIENT
        constexpr int jacobianHeight = 16;
        outputs[0] = post.cl;
        outputs[1] = post.cd;
        for (int i=0;i<7;++i){
            outputs[2+i] = topsurf[i];
            outputs[2+7+i] = botsurf[i];
        }
        #elif DO_SOUND
        constexpr int jacobianHeight = 1;
        outputs = OASPL;
        #else
        constexpr int jacobianHeight = 2;
        outputs[0] = post.cl;
        outputs[1] = post.cd;
        #endif

        #if DO_SOUND
        tape.registerOutput(outputs);
        tape.setPassive();
        outputs.gradient() = 1.0 ;
        #else
        for (int i=0;i<jacobianHeight;++i){
            tape.registerOutput(outputs[i]);
        }
        tape.setPassive();
        for (int i=0;i<jacobianHeight;++i){outputs[i].gradient()[i] = 1.0 ;}
        #endif
        
        tape.evaluate();

        codi::Jacobian<double> jacobian(jacobianHeight,Nin);
        codi::Jacobian<double> jacobianAlpha(jacobianHeight,1);
        codi::Jacobian<double> jacobianRe(jacobianHeight,1);
        codi::Jacobian<double> jacobianMa(jacobianHeight,1);
        
        #if DO_SOUND
        std::vector<std::vector<double>> allGradients ;
        for (int i = 0; i < Nin; ++i) {    
            jacobian(0,i) = inYcoords[i].getGradient();
        }
        jacobianAlpha(0,0) = alphad.getGradient();
        jacobianRe(0,0)    = Re.getGradient();
        jacobianMa(0,0)    = Ma.getGradient();
        #else
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
        #endif
        


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
            j["d " + outputNames[i] + " / d Re"] = allGradientsRe[i];
            j["d " + outputNames[i] + " / d Ma"] = allGradientsMa[i];
        }
        
        #if DO_SOUND
        std::ofstream outFile("OASPL_gradients.json");
        #else
        std::ofstream outFile("ad_gradients.json");
        #endif
        outFile << j.dump(4);  // pretty-print with 4-space indentation
        outFile.close();
        
        tape.reset();

    #endif

    return converged;
    //}
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
    Real customUinf = j["Uinf"].get<double>();
    int useCustUinf = j["custUinf"].get<int>();
    const Real X = j["X"].get<double>();
    const Real Y = j["Y"].get<double>();
    const Real Z = j["Z"].get<double>();
    const Real S = j["S"].get<double>();
    const Real Ncrit = j["ncrit"].get<double>();
    const int doCps = j["returnData"].get<int>();
    const Real Ufac = j["Ufac"].get<double>();
    const Real TEfac = j["TEfac"].get<double>();

    // forcing transition
    const bool force = j["forcetrans"].get<int>();
    const Real topTransPos = j["toptrans"].get<double>();
    const Real botTransPos = j["bottrans"].get<double>();
    
    const bool custNcrits =  j["custncrits"].get<int>();
    const Real topNcrit = j["topncrit"].get<double>();
    const Real botNcrit = j["botncrit"].get<double>();


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
    

    bool converged = runCode(doRestart,doXfoilStart,doGetPoints,targetAlphaDeg,Re,Ma,inXcoords,inYcoords,initStates,initTurb,sampleTE,X,Y,Z,S,customUinf,useCustUinf,doCps,Ncrit,Ufac,TEfac,topTransPos,botTransPos,force,custNcrits,topNcrit,botNcrit);
    
    std::cout << "converged: " << converged << std::endl ;
    return converged;
};