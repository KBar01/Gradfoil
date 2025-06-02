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

def use_xfoil(xcoords,ycoords,alphaDeg,Re=1e6,Ma=0.0,sampleTE=0.95,xfoilPath=None):
    
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords": xcoords,
        "ycoords": ycoords,
        "alpha_degrees": alphaDeg,
        "Re": Re,
        "Ma": Ma,
        "restart": 0,
        "xfoilstart":0,
        "xfoilgetpoints":0,
        "sampleTE": sampleTE
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)
    
    
    completed = xfoil_start_run(alphaDeg,Re,Ma,xcoords,ycoords,sampleTE,EXEC_FWD,xfoilPath)
    return completed

def standard_run(xcoords,ycoords,alphaDeg,Re=1e6,Ma=0.0,sampleTE=0.95,xfoilPath=None):
    
    
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords": xcoords,
        "ycoords": ycoords,
        "alpha_degrees": alphaDeg,
        "Re": Re,
        "Ma": Ma,
        "restart": 0,
        "xfoilstart":0,
        "xfoilgetpoints":0,
        "sampleTE": sampleTE
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting
    initResult = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode

    completed = 0

    stepsize = 1.0
    if (not initConvergence):


        # Starting back stepping to lower AoA, searching for converged solution 
        max_back_steps = 10
        min_alpha = -1.0
        back_converged = False
        tempalf = np.round(alphaDeg,decimals=1) - stepsize ; 
        smallStep = 0.5
        for i in range(max_back_steps):
            
            if tempalf < 2 :
                tempalf = np.round(tempalf,decimals=1)
                smallStep = 0.1
            
            if tempalf < min_alpha:
                print("Minimum AoA reached, cannot backstep further.")
                break

            #print(f"Trying backstep to: {tempalf:.17f} degrees")
            
            with open(in_json_path, "r") as f:
                data = json.load(f)
            data["alpha_degrees"] = tempalf
            data["restart"] = 0
            with open(in_json_path, "w") as f:
                json.dump(data, f)
            
            backstepResult = subprocess.run([EXEC_FWD], cwd=os.getcwd(), capture_output=True, text=True)
            
            if backstepResult.returncode == 1:
                print(f"Backstep converged at {tempalf}")
                back_converged = True
                break
            
            
            tempalf -= smallStep

        # At this point have stepped beckward to a converged solution for given aerofoil
        # now need to step forward to reach the initial target alpha

        if not back_converged:
            return False
        
        
        stepsize = 0.5
        fwdalf = tempalf + stepsize
        attemptCount = 0
        
        overallCount = 0
        while not completed:
            
            # attempt fwd step using prev solution written to restart file
            print(f"Trying forward step to: {fwdalf} degrees")


            with open(in_json_path, "r") as f:
                    data = json.load(f)

            data["restart"] = 1 # using restart fle
            data["alpha_degrees"] = fwdalf  # trying different alfa
            
            with open(in_json_path, "w") as f:
                json.dump(data, f, indent=4)  # indent is optional, just for readability

            result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
            fwdStepConverged = result.returncode
            
            
            if fwdStepConverged:
                if abs(fwdalf - alphaDeg) < 1e-3:
                    completed = True
                    break
                else:
                    fwdalf = min(fwdalf + stepsize, alphaDeg)
                    attemptCount = 0  # reset on success
            else:
                attemptCount += 1
                if attemptCount > 6:
                    print("Forward stepping failed after multiple attempts")
                    break
                fwdalf -= stepsize / (2 ** attemptCount)
        
            if overallCount > 15:
                break
            overallCount +=1
        
        
        if (not completed) and (xfoilPath != None):

            # c++ code not converging, use xfoil to give an initial starting point
            completed = xfoil_start_run(alphaDeg,Re,Ma,xcoords,ycoords,EXEC_FWD,xfoilPath)
        
        return completed
    

    return (initConvergence or completed)


def fwd_run(xcoords,ycoords,alphaDeg,Re=1e6,Ma=0.0,sampleTE=0.95,xfoilPath=None,xfoilStart=0):
    
    success = False
    if xfoilStart == 0:
        success = standard_run(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,xfoilPath)
    else:
        success = use_xfoil(xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,xfoilPath)
    return success


def grad_run():

    # Run the AD version of the code, using known solution from fwd run
    result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)




