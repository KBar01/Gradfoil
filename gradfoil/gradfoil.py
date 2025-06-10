import json
import subprocess
import numpy as np
import os
import tempfile
from importlib.resources import files
from .xfoilExtract import xfoil_start_run

# Dynamically locate the installed executable path (in gradfoil/bin/)
BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
EXEC_FWD = os.path.join(BIN_DIR, "CFoil_fwd")
EXEC_AD = os.path.join(BIN_DIR, "CFoil_AD")
EXEC_NOISE = os.path.join(BIN_DIR, "CFoil_Noise")
def use_xfoil(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,X,Y,Z,S,xfoilPath,Uinf,custUinf,trackCLs,ncrit):
    
    completed = xfoil_start_run(alphaDeg,Re,Ma,xcoords,ycoords,sampleTE,X,Y,Z,S,EXEC_FWD,xfoilPath,Uinf,custUinf,trackCLs,ncrit)
    return completed

def standard_run(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,X,Y,Z,S,xfoilPath,Uinf,custUinf,trackCLs,returnFoilCps,ncrit,Ufac,TEfac):
    
    
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords":       xcoords,
        "ycoords":       ycoords,
        "alpha_degrees": alphaDeg,
        "Re":            Re,
        "Ma":            Ma,
        "restart":       0,
        "xfoilstart":    0,
        "xfoilgetpoints":0,
        "sampleTE":      sampleTE,
        "X":             X,
        "Y":             Y,
        "Z":             Z,
        "S":             S,
        "Uinf":          Uinf,
        "custUinf":      custUinf,
        "returnData":    returnFoilCps,
        "ncrit":         ncrit,
        "Ufac":          Ufac,
        "TEfac":         TEfac
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting
    initResult = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode

    if initConvergence:
        return True

    print("Initial run failed. Starting backstepping ...")
    
    max_back_steps = 10
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

        # Attempt run
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
    max_attempts = 10

    while (not completed) and (overallCount <= max_attempts):
        
        print(f"Trying forward step to: {fwdalf:.2f}")

        with open(in_json_path, "r") as f:
            data = json.load(f)

        data["restart"] = 1
        data["alpha_degrees"] = fwdalf
        with open(in_json_path, "w") as f:
            json.dump(data, f, indent=4)

        result = subprocess.run([EXEC_FWD], cwd=os.getcwd(), capture_output=True, text=True)
        converged = result.returncode == 1

        if converged:
            if abs(fwdalf - alphaDeg) < 1e-3:
                completed = True
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

    # Optional fallback to XFOIL
    if (not completed) and (xfoilPath is not None):
        print("Falling back to XFOIL for initialization...")
        completed = xfoil_start_run(alphaDeg,Re,Ma,xcoords,ycoords,sampleTE,X,Y,Z,S,EXEC_FWD,xfoilPath,Uinf,custUinf,trackCLs,ncrit)

    return completed
    



def fwd_run(xcoords,ycoords,alphaDeg,Re=1e6,Ma=0.0,sampleTE=0.95,observerX=0.0,observerY=0.0,observerZ=1.2,span=0.5,xfoilPath=None,xfoilStart=0,Uinf=1,custUinf=0,trackCLs=0,returnFoilCps=0,ncrit=9.0,Ufac=2.5,TEfac=0.06):
    
    success = False
    if xfoilStart == 0:
        success = standard_run(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,observerX,observerY,observerZ,span,xfoilPath,Uinf,custUinf,trackCLs,returnFoilCps,ncrit,Ufac,TEfac)

        if success:
            return success, Ufac,TEfac
    
        if not success:
            for uf, tef in [(Ufac-0.5, TEfac+0.01), (Ufac+0.5, TEfac-0.02), (Ufac+0.1, TEfac+0.05)]:
                
                print('trying different panel distribution')
                success = standard_run(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,observerX,observerY,observerZ,span,xfoilPath,Uinf,custUinf,trackCLs,returnFoilCps,ncrit,uf,tef)
                if success:
                    return success, uf,tef
        
        return success,Ufac,TEfac

    else:
        success =    use_xfoil(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,observerX,observerY,observerZ,span,xfoilPath,Uinf,custUinf,trackCLs,ncrit)
        return success,0.0,0.0


def grad_run(doSound=0):
    # Run the AD version of the code, using known solution from fwd run
    if doSound:
        result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)
        result = subprocess.run([EXEC_NOISE],cwd=os.getcwd(), capture_output=True, text=True)

    else:
        result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)




