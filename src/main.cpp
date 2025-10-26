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
    bool fwdDoubleRestart,
    bool fwdCodiRestart,
    bool ADrestart,
    bool doLateADRestart,
    const Real ADRestartRnorm,
    Real alphad,
    Real Re, 
    Real Ma,
    Real rhoInf,
    Real kinViscInf,
    const Real (&inXcoords)[Nin], 
    Real (&inYcoords)[Nin],
    const Real sampleTE,
    const Real X,
    const Real Y,
    const Real Z,
    const Real S,
    Real Uinf,
    const int doCps,
    const Real nCrit,
    const Real Ufac, 
    const Real TEfac,
    const Real &topTransPos,
    const Real &botTransPos,
    const bool force,
    const std::string model){

    
    Real outputs[3] ;

    // ------------- Doing Adjoint, register relevant input to track gradients --------------------------
    #ifdef AD_VERSION

        Tape& tape = Real::getTape();
        tape.setActive();

        for (int i = 0; i < Nin; ++i) {
            tape.registerInput(inYcoords[i]);
        }
        tape.registerInput(alphad);

    #elif FWD_CODI_VERSION

        //Tape& tape = Real::getTape();
        //tape.setActive();
        //tape.registerInput(X);
        //tape.setPassive();

    #endif

    Real alpha = (alphad/180)*M_PI;
    
    Real minX = 0.5, maxX = 0.01;
    for (int i=0;i<Nin;++i){
        Real newMin = std::min(minX,inXcoords[i]);
        Real newMax = std::max(maxX,inXcoords[i]);
        minX = newMin ;
        maxX = newMax ;
    }
    Real chordScale = maxX - minX ;
    
    Oper oper(alpha,Re,Ma);
    oper.rho = rhoInf;

    
    if (Uinf <= 0.0){   // indicates no custom Uinf given
        Uinf = (Re*kinViscInf)/(chordScale) ;
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

    Foil foil(flattenedCoords);
    Isol isol;
    Param param;
    param.ncrit = nCrit;
    Wake wake;
    Vsol vsol;
    
    static Glob glob;

    glob.doADrestartExtract = doLateADRestart;
    glob.ADrestartRnorm = ADRestartRnorm ;
    
    build_gamma_codi(isol,foil,oper);
    init_thermo(oper,param,geom);
    build_wake(foil,geom,oper,isol,wake);
    stagpoint_find(isol,foil,wake);
    identify_surfaces(isol,vsol);
    set_wake_gap(foil,isol,vsol);
    calc_ue_m(foil,wake,isol,vsol);
    rebuild_ue_m(foil,wake,isol,vsol,false);
    
    #ifdef AD_VERSION
    
        if (ADrestart){

            std::ifstream prevfile("prevRestart.json");
            if (!prevfile) {
                std::cerr << "Failed to open restart.json\n";
                return 1;
            }

            // Parse the JSON
            json j;
            prevfile >> j;

            for (int i = 0; i < RVdimension; ++i) {
            double val = j["states"][i].get<double>();
            glob.U[i] = val;  // assigns numeric value to RealReverse
            }

            for (int i = 0; i < (Ncoords + Nwake); ++i) {
            double val = j["turb"][i].get<double>();
            vsol.turb[i] = val;
            }

            if (force){
                tdata.isForced[0]  = 1;
                tdata.isForced[1]  = 1;
            }
        }
        else {
            init_boundary_layer(oper,foil,param,isol,vsol,glob,tdata,force);
        }
        
    #elif FWD_DOUBLE_VERSION

        // fwd_double run, dealing with restart from different alpha or not
        if (fwdDoubleRestart){

            std::ifstream prevfile("prevRestart.json");
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
        else {
            init_boundary_layer(oper,foil,param,isol,vsol,glob,tdata,force);
        }

    #else  // doing FWD_CODI_VERSION

        if (fwdCodiRestart){

            std::ifstream prevfile("prevRestart.json");
            if (!prevfile) {
                std::cerr << "Failed to open restart.json\n";
                return 1;
            }

            // Parse the JSON
            json j;
            prevfile >> j;

            for (int i = 0; i < RVdimension; ++i) {
            double val = j["states"][i].get<double>();
            glob.U[i] = val;  // assigns numeric value to RealReverse
            }

            for (int i = 0; i < (Ncoords + Nwake); ++i) {
            bool val = j["turb"][i].get<bool>();
            vsol.turb[i] = val;
            }

            if (force){
                tdata.isForced[0]  = 1;
                tdata.isForced[1]  = 1;
            }
        }
        else {
            init_boundary_layer(oper,foil,param,isol,vsol,glob,tdata,force);
        }
    #endif

    stagpoint_move(isol,glob,foil,wake,vsol);
    bool converged = solve_coupled(oper,foil,wake,param,vsol,isol,glob,tdata,force);
    Post post;
    calc_force(oper,geom,param,isol,foil,glob,post);

    #ifndef AD_VERSION
        Real tauWall[Ncoords];
        if (doCps){
            Real cf_U[4]={0};
            for (int i=0;i<Ncoords;++i){
                tauWall[i] = get_cf(glob.U[colMajorIndex(0,i,4)],
                                    glob.U[colMajorIndex(1,i,4)],
                                    glob.U[colMajorIndex(2,i,4)],
                                    glob.U[colMajorIndex(3,i,4)],
                                    vsol.turb[i],
                                    false,
                                    param,
                                    cf_U
                );
                
                tauWall[i] *=  (oper.rho*(glob.U[colMajorIndex(3,i,4)]*Uinf * glob.U[colMajorIndex(3,i,4)]*Uinf))/2;
            }
        }
    #endif


    Real topsurf[7],botsurf[7];
    Real xcoords[Ncoords]={0};
    Real ycoords[Ncoords]={0};

    for (int i=0;i<Ncoords;++i){
        xcoords[i] = flattenedCoords[colMajorIndex(0,i,2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1,i,2)];
    }

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,vsol,param,topsurf,botsurf,Uinf,geom,(sampleTE*geom.chord));

    Real OASPL = calc_OASPL(botsurf,topsurf,oper,geom,Uinf,X,Y,Z,S,kinViscInf,oper.rho,doCps,model);
    
    std::vector<std::string> outputNames = {"CL", "CD", "OASPL"};
    

    # ifdef FWD_CODI_VERSION
    
        // check OASPL validity
        if (std::isnan(OASPL) || std::isinf(OASPL)) {
            converged = false;
        }
        
        if (converged){

            json out;
            out["aerofoilChord"] = chordScale.getValue();
            out["freestreamVelocity"] = Uinf.getValue();
            out["CL"]  = post.cl.getValue();
            out["CD"]  = post.cd.getValue();
            out["conv"] = 1;
            out["OASPL"] = OASPL.getValue();
            
            
            if (doCps){
                
                json restart;
                std::vector<double> states_d(RVdimension);
                // Extract numeric values from CoDiPack types
                for (size_t i = 0; i < RVdimension; ++i){
                    states_d[i] = glob.U[i].getValue();
                } 
                restart["states"] = states_d;
                restart["turb"]   = vsol.turb;
                std::ofstream restartFile("restart.json");
                restartFile << restart.dump(4);  // pretty print with 4 spaces indentation
                restartFile.close();

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

                double inner[2*Ncoords] ;
                double cps[2*Ncoords];
                for (int i=0;i<(2*Ncoords);++i){
                    inner[i] = foil.x[i].getValue() ;
                    cps[i] = post.cp[i].getValue();
                }
                out["innerFoil"] = inner;
                out["Cp"] = cps;
                
                
                out["stagnation"] = isol.stagIndex;
                out["topTransX"]  = topTransX.getValue();
                out["botTransX"]  = botTransX.getValue();
                
                
                //out["tauWall"] = tauWall;
                std::vector<std::string> BLoutputNames = {"CL", "CD",
                    "thetaUpper", "deltaStarUpper", "tauMaxUpper","edgeVelocityUpper", "dpdxUpper", "tauWallUpper", "delta99Upper",
                    "thetaLower", "deltaStarLower", "tauMaxLower","edgeVelocityLower", "dpdxLower", "tauWallLower", "delta99Lower"
                };
                out[BLoutputNames[2]] = topsurf[0].getValue();
                out[BLoutputNames[3]] = topsurf[1].getValue();
                out[BLoutputNames[4]] = topsurf[2].getValue();
                out[BLoutputNames[5]] = topsurf[3].getValue();
                out[BLoutputNames[6]] = topsurf[4].getValue();
                out[BLoutputNames[7]] = topsurf[5].getValue();
                out[BLoutputNames[8]] = topsurf[6].getValue();

                out[BLoutputNames[9]] = botsurf[0].getValue();
                out[BLoutputNames[10]] = botsurf[1].getValue();
                out[BLoutputNames[11]] = botsurf[2].getValue();
                out[BLoutputNames[12]] = botsurf[3].getValue();
                out[BLoutputNames[13]] = botsurf[4].getValue();
                out[BLoutputNames[14]] = botsurf[5].getValue();
                out[BLoutputNames[15]] = botsurf[6].getValue();
            }
            
            std::ofstream outFile("out.json");
            outFile << out.dump(4);  // pretty print with 4 spaces indentation
            outFile.close();
            
            
        }
        else{

            json out;
            out["conv"] = 0;
            std::ofstream outFile("out.json");
            outFile << out.dump(4);  // pretty print with 4 spaces indentation
            outFile.close();

        }

    #elif FWD_DOUBLE_VERSION
        // check OASPL validity
        if (std::isnan(OASPL)   // caught NaNs
            || std::isinf(OASPL)) // caught infs     
        {
            converged = false;
        }
    
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

            out["conv"] = 1;
            out["OASPL"] = OASPL;

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
                
                
                out["tauWall"] = tauWall;
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
   
    #ifdef AD_VERSION

        
        constexpr int jacobianHeight = 3;
        outputs[0] = post.cl;
        outputs[1] = post.cd;
        outputs[2] = OASPL;
            

        for (int i=0;i<jacobianHeight;++i){
            tape.registerOutput(outputs[i]);
        }
        tape.setPassive();
        for (int i=0;i<jacobianHeight;++i){outputs[i].gradient()[i] = 1.0 ;}
        
        tape.evaluate();

        codi::Jacobian<double> jacobian(jacobianHeight,Nin);
        codi::Jacobian<double> jacobianAlpha(jacobianHeight,1);

        
        std::vector<std::vector<double>> allGradients ;
        for (int i = 0; i < Nin; ++i) {   
            for (int n=0;n<jacobianHeight;++n){
                jacobian(n,i) = inYcoords[i].getGradient()[n];
            }
        }
        for (int n=0;n<jacobianHeight;++n){
            jacobianAlpha(n,0) = alphad.getGradient()[n];
        }
        
        // Fill allGradients
        for (int out = 0; out < jacobianHeight; ++out) {
            std::vector<double> grad;
            for (int i = 0; i < Nin; ++i) {
                grad.push_back(jacobian(out, i));
            }
            allGradients.push_back(grad);
        }
        // Create JSON
        json j;
        for (int i = 0; i < allGradients.size(); ++i) {
            j["d " + outputNames[i] + " / d ycoords"] = allGradients[i];
        }
        
        std::vector<double> allGradientsAlf ;

        for (int out = 0; out<jacobianHeight; ++out) {
            allGradientsAlf.push_back(jacobianAlpha(out, 0));
        }
        
        for (int i = 0; i < allGradientsAlf.size(); ++i) {
            j["d " + outputNames[i] + " / d alpha"] = allGradientsAlf[i];
        }
        
        std::ofstream outFile("ad_gradients.json");
        outFile << j.dump(4);  // pretty-print with 4-space indentation
        outFile.close();
        
        json js;

        js["CL"] = post.cl.getValue();
        js["CD"] = post.cd.getValue();

        std::ofstream newoutFile("ad_outs.json");
        newoutFile << js.dump(4);
        newoutFile.close();
        
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

    Real inXcoords[Nin]={0}, inYcoords[Nin]={0};

    for (int i = 0; i < Nin; ++i) {
        inXcoords[i] = j["xcoords"][i];       // X[0][i] -> x[i] (x-coordinates)
        inYcoords[i] = j["ycoords"][i];   // X[1][i] -> x[N + i] (y-coordinates)
    }
   

    // Read input variables
    Real targetAlphaDeg = j["alpha_degrees"].get<double>();
    Real Re = j["Re"].get<double>();
    Real Ma = j["Ma"].get<double>();
    Real rhoInf = j["rho"].get<double>();
    Real nuInf = j["nu"].get<double>();

    int doRestart = j["restart"].get<int>();
    int fwdCodiRestart = j["fwdCodiRestart"].get<int>();
    int doADRestart = j["ADrestart"].get<int>();
    int doLaterestart = j["lateRestart"].get<int>();
    const Real lateRestartNorm = j["lateRestartNorm"].get<double>();


    Real sampleTE = j["sampleTE"].get<double>();
    Real customUinf = j["Uinf"].get<double>();
    const Real X = j["X"].get<double>();
    const Real Y = j["Y"].get<double>();
    const Real Z = j["Z"].get<double>();
    const Real S = j["S"].get<double>();
    const Real Ncrit = j["ncrit"].get<double>();
    const int doCps = j["returnData"].get<int>();
    const Real Ufac = j["Ufac"].get<double>();
    const Real TEfac = j["TEfac"].get<double>();

    // forcing transition variables
    const bool force = j["forcetrans"].get<int>();
    const Real topTransPos = j["toptrans"].get<double>();
    const Real botTransPos = j["bottrans"].get<double>();
    
    const std::string model = j["model"].get<std::string>();

   
    bool converged = runCode(doRestart,fwdCodiRestart,doADRestart,doLaterestart,lateRestartNorm,
                            targetAlphaDeg,Re,Ma,rhoInf,nuInf,inXcoords,inYcoords,
                            sampleTE,X,Y,Z,S,customUinf,
                            doCps,Ncrit,Ufac,TEfac,
                            topTransPos,botTransPos,force,
                            model);
    
    return converged;
    
};