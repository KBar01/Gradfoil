import json
import subprocess
import numpy as np
import os
#import tempfile
#from importlib.resources import files
from .xfoilExtract import xfoil_start_run

# Dynamically locate the installed executable path (in gradfoil/bin/)
BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
EXEC_FWD = os.path.join(BIN_DIR, "CFoil_fwd_double")
EXEC_FWD_codi = os.path.join(BIN_DIR, "CFoil_fwd_codi")
EXEC_AD = os.path.join(BIN_DIR, "CFoil_AD")


def standard_run(xcoords,ycoords,Re,alphaDeg,Ma,sampleTE,X,Y,Z,S,model,rho,nu,ncrit,custUinf,returnFoilCps,Ufac,TEfac,toptrans,bottrans,force,lateRestart,lateRestartNorm,chord):
    

    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords":       xcoords,
        "ycoords":       ycoords,
        "alpha_degrees": alphaDeg,
        "Re":            Re,
        "Ma":            Ma,
        "rho":           rho,
        "nu":            nu,
        "restart":        0,
        "fwdCodiRestart": 0,
        "ADrestart":      0,
        "lateRestart":    lateRestart,
        "lateRestartNorm": lateRestartNorm,
        "sampleTE":      sampleTE,
        "X":             X,
        "Y":             Y,
        "Z":             Z,
        "S":             S,
        "Uinf":      custUinf,
        "returnData":    returnFoilCps,
        "ncrit":         ncrit,
        "Ufac":          Ufac,
        "TEfac":         TEfac,
        "toptrans":      toptrans,
        "bottrans":      bottrans,
        "forcetrans":    force,
        "model":  model,
        "chord": chord,
        "WPSonly":0
    }

    # Write JSON input file
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting, use codi version to ensure output match to AD version of code
    initResult = subprocess.run([EXEC_FWD_codi],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode
    if initConvergence==1:
        return True

    print("Initial run failed. Starting backstepping ...")
    
    max_back_steps = 5
    stepsize = 1.0
    small_step = 0.5
    back_converged = False
    completed = False
    
    # Determine stepping direction based on sign of alphaDeg
    if alphaDeg >= 0:
        step_direction = -1
    else:
        step_direction = 1
    
    # take a step of 0.5degrees, set limit of backstep to start + 5.0
    tempalf = np.round(alphaDeg, decimals=1) + (step_direction * stepsize)
    min_alpha = alphaDeg + (step_direction * 5.0)
    
    for i in range(max_back_steps):
        
        if abs(tempalf) < 2.0:
            small_step = 0.1

        # Check if minimum/maximum alpha reached
        if (step_direction < 0 and tempalf < min_alpha) or (step_direction > 0 and tempalf > min_alpha):
            print("Minimum backstep AoA reached. Cannot continue.")
            break

        # Modify input JSON
        with open(in_json_path, "r") as f:
            data = json.load(f)

        data["alpha_degrees"] = tempalf
        data["restart"] = 0  # fresh run
        with open(in_json_path, "w") as f:
            json.dump(data, f)

        # Attempt run, use fwd_double version to be quicker 
        result = subprocess.run([EXEC_FWD], cwd=os.getcwd(), capture_output=True, text=True)
        if result.returncode == 1:
            print(f"Backstep converged at {tempalf}")
            back_converged = True
            break

        tempalf += (step_direction * small_step)

    if not back_converged:
        print("Backstepping failed. No converged base solution.")
        return False
    
    
    # Step forward toward original alphaDeg using restart
    print("Starting forward stepping...")
    stepsize = 0.5
    fwdalf = tempalf - (step_direction * stepsize)
    attemptCount = 0
    overallCount = 0
    max_attempts = 6

    while (not completed) and (overallCount <= max_attempts):
        
        print(f"Trying forward step to: {fwdalf:.2f}")

        with open(in_json_path, "r") as f:
            data = json.load(f)

        data["restart"] = 1
        data["alpha_degrees"] = fwdalf
        with open(in_json_path, "w") as f:
            json.dump(data, f, indent=4)

        # Rename restart.json to prevRestart.json before running the executable
        restart_path = os.path.join(os.getcwd(), "restart.json")
        prev_restart_path = os.path.join(os.getcwd(), "prevRestart.json")

        if os.path.exists(restart_path):
            # If a previous prevRestart.json exists, remove it to avoid clutter
            if os.path.exists(prev_restart_path):
                os.remove(prev_restart_path)
            os.rename(restart_path, prev_restart_path)

        result = subprocess.run([EXEC_FWD], cwd=os.getcwd(), capture_output=True, text=True)
        converged = result.returncode == 1

        if converged:
            if abs(fwdalf - alphaDeg) < 1e-3:
                
                with open(in_json_path, "r") as f:
                    data = json.load(f)
                data["fwdCodiRestart"] = 1
                data["restart"] = 0
                with open(in_json_path, "w") as f:
                    json.dump(data, f, indent=4)
                
                result = subprocess.run([EXEC_FWD_codi], cwd=os.getcwd(), capture_output=True, text=True)

                with open(in_json_path, "r") as f:
                    data = json.load(f)
                data["restart"] = 0
                data["fwdCodiRestart"] = 0
                data["ADrestart"] = 1
                with open(in_json_path, "w") as f:
                    json.dump(data, f, indent=4)
                    converged = result.returncode == 1
                    completed = converged
                
                break
            else:
                
                nextStep = fwdalf - (step_direction*stepsize)
                diff = alphaDeg - fwdalf

                if abs(diff) < abs(nextStep):
                    fwdalf = alphaDeg
                else:
                    fwdalf = nextStep

                attemptCount = 0
        else:
            attemptCount += 1
            if attemptCount > 6:
                print("Forward stepping failed repeatedly.")
                break

            fwdalf += step_direction * (stepsize / (2 ** attemptCount))

        overallCount += 1

    return completed
    



def fwd_run(xcoords,ycoords,alphaDeg,Re=1e6,Ma=0.0,
            sampleTE=0.95,observerX=0.0,observerY=0.0,observerZ=1.2,span=0.5, model="kam",
            rho=1.225, nu=1.789e-5, Uinf=0,ncrit=9.0,
            Ufac=1.0,TEfac=0.09, repanel=0,
            toptrans=0.5,bottrans=0.5,forcetrans=0,
            returnFoilCps=0,lateRestart=0,lateRestartNorm=1e-5,chord=1.0):
    
    
    #xcoords,ycoords,Re,alphaDeg,Ma,sampleTE,X,Y,Z,S,model,rho,nu,ncrit,custUinf,returnFoilCps,Ufac,TEfac,toptrans,bottrans,force,lateRestart,lateRestartNorm
    if repanel:
        success = standard_run(xcoords,ycoords,Re,alphaDeg,Ma,
                               sampleTE,observerX,observerY,observerZ,span,model,
                               rho,nu,ncrit,Uinf,returnFoilCps,Ufac,TEfac,toptrans,bottrans,forcetrans,lateRestart,lateRestartNorm,chord)

        if success:
            return success
        else:
            count = 1 
            for uf, tef in [(1.8,0.1), (2.1,0.09), (2.6,0.09), (1.0,0.09), (1.0,1.1), (1.5,0.09)]:
                
                print('trying different panel distribution ('+str(count)+'/6)')
                success = standard_run(xcoords,ycoords,Re,alphaDeg,Ma,
                               sampleTE,observerX,observerY,observerZ,span,model,
                               rho,nu,ncrit,Uinf,returnFoilCps,uf,tef,toptrans,bottrans,forcetrans,lateRestart,lateRestartNorm,chord)
                if success:
                    break
                count +=1
        
            return success
    
    else:
        success = standard_run(xcoords,ycoords,Re,alphaDeg,Ma,
                               sampleTE,observerX,observerY,observerZ,span,model,
                               rho,nu,ncrit,Uinf,returnFoilCps,Ufac,TEfac,toptrans,bottrans,forcetrans,lateRestart,lateRestartNorm,chord)
        return success
    


def grad_run():
    # Run the AD version of the code, using known solution from fwd run
    result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)

def WPS_run(Re,X,Y,Z,S,model,chord,rho,nu,topBLstates,botBLstates):
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "Re":            Re,
        "rho":           rho,
        "nu":            nu,
        "X":             X,
        "Y":             Y,
        "Z":             Z,
        "S":             S,
        "model":  model,
        "chord": chord,
        "WPSonly": 1,
        
        # -------- TOP boundary layer states --------
        # order: dstar, theta, delta, tauw, taumax, ue, dpdx
        "topdstar":   topBLstates[0],
        "toptheta":   topBLstates[1],
        "topdelta":   topBLstates[2],
        "toptauw":    topBLstates[3],
        "toptaumax":  topBLstates[4],
        "topue":      topBLstates[5],
        "topdpdx":    topBLstates[6],

        # -------- BOTTOM boundary layer states --------
        "botdstar":   botBLstates[0],
        "bottheta":   botBLstates[1],
        "botdelta":   botBLstates[2],
        "bottauw":    botBLstates[3],
        "bottaumax":  botBLstates[4],
        "botue":      botBLstates[5],
        "botdpdx":    botBLstates[6],
    }

    # Write JSON input file
    with open(in_json_path, "w") as f:
        json.dump(data, f,indent=4)

    # Run the executable for first time, no restarting, use codi version to ensure output match to AD version of code
    initResult = subprocess.run([EXEC_FWD_codi],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode
    if initConvergence==1:
        return True
    
    





