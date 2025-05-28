import json
import subprocess
import numpy as np
import os
import tempfile
from importlib.resources import files
from xfoilRun import run_xfoil_get_BL_states

# Dynamically locate the installed executable path (in gradfoil/bin/)
BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
EXEC_FWD = os.path.join(BIN_DIR, "CFoil_fwd")
EXEC_AD = os.path.join(BIN_DIR, "CFoil_AD")


def standard_run(alphaDeg,Re,Ma,xcoords,ycoords):
    
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords": xcoords,
        "ycoords": ycoords,
        "alpha_degrees": alphaDeg,
        "Re": Re,
        "Ma": Ma,
        "restart": 0,
        "xfoilStart":0
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting
    initResult = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode

    completed = 0

    stepsize = 0.5
    if (not initConvergence):

        max_back_steps = 50
        min_alpha = 0.0
        back_converged = False
        tempalf = np.round(alphaDeg,decimals=1) - stepsize ; 
        smallStep = 0.25
        for i in range(max_back_steps):
            
            if tempalf < 2.05 :
                tempalf = np.round(tempalf,decimals=1)
                smallStep = 0.1
            
            if tempalf < min_alpha:
                print("Minimum AoA reached, cannot backstep further.")
                break

            print(f"Trying backstep to: {tempalf:.17f} degrees")
            
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
                else:
                    fwdalf = min(fwdalf + stepsize, alphaDeg)
                    attemptCount = 0  # reset on success
            else:
                attemptCount += 1
                if attemptCount > 6:
                    print("Forward stepping failed after multiple attempts")
                    break
                fwdalf -= stepsize / (2 ** attemptCount)
        
            if overallCount > 20:
                break
            overallCount +=1
        
        
        return completed
    

    return (initConvergence or completed)


def grad_run():

    # Run the AD version of the code, using known solution from fwd run
    result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)




def xfoil_start_run(alphaDeg,Re,Ma,xcoords,ycoords):
    

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
        "xfoilgetpoints":1
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable to get inner foil coords
    initResult = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode

    foil_json_path = os.path.join(cwd, "innerfoilNodes.json")
    with open(foil_json_path, "r") as f:
        data = json.load(f)
    
    nodes = data["points"]
    nodes.reshape((2, 200), order='F')

    # Run xfoil with the foil for converged states
    innerFoil = nodes.T
    xfoilStates,xfoilturb = run_xfoil_get_BL_states(innerFoil,alphaDeg,Re,Ma)
    # save them to json
    xfoil_json_path = os.path.join(cwd, "xfoilstates.json")
    data = {
        "states":xfoilStates.flatten(order='F'),
        "turb": xfoilturb
    }
    with open(xfoil_json_path, "w") as f:
        json.dump(data, f)

    with open(in_json_path, "r") as f:
        data = json.load(f)
    data["xfoilgetpoints"] = 0 
    data["xfoilstart"] = 1      
    with open(in_json_path, "w") as f:
        json.dump(data, f, indent=4)
    
    Result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    converged = Result.returncode

    return converged
