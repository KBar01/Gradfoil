import json
import subprocess
import numpy as np
import os
import tempfile
from importlib.resources import files

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
        "restart": 0
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting
    result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = result.returncode

    stepsize = 0.5
    if (not initConvergence):
        tempalf = alphaDeg - stepsize ; 
        for i in range(50):
            
            with open(in_json_path, "r") as f:
                data = json.load(f)

            # Modify only the desired fields
            #data["restart"] = 1  # or 0
            data["alpha_degrees"] = tempalf  # assign your new value here

            # Write back the updated JSON
            with open(in_json_path, "w") as f:
                json.dump(data, f) 


            result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
            backStepConverged = result.returncode

            if (backStepConverged):
                print("stepped back to: " + str(tempalf))
                break

            tempalf -= 0.5

        # At this point have stepped beckward to a converged solution for given aerofoil
        # now need to step forward to reach the initial target alpha

        fwdalf = tempalf + stepsize
        attemptCount = 0
        completed = 0
        while ((not completed)):
            
            # attempt fwd step using prev solution written to restart file

            print("solving : " +str(fwdalf))

            with open(in_json_path, "r") as f:
                    data = json.load(f)

            data["restart"] = 1 # using restart fle
            data["alpha_degrees"] = fwdalf  # trying different alfa
            
            with open(in_json_path, "w") as f:
                json.dump(data, f, indent=4)  # indent is optional, just for readability

            result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
            fwdStepConverged = result.returncode
            if (not fwdStepConverged):
                
                if attemptCount > 6 :
                    break
                attemptCount += 1
                fwdalf -= (stepsize / (2**attemptCount)) # take smaller and smaller steps
                continue
            
            if (alphaDeg - fwdalf > (stepsize + 0.01)):
                fwdalf += stepsize              # step fwd 0.25
            elif (fwdalf == alphaDeg):
                break                       # just solved target alf, break loop
            else:
                fwdalf = alphaDeg # step fwd to target alf
            
            if (alphaDeg - fwdalf > stepsize):
                attemptCount = 0;  # reset attempt as moving forward with new solution at higher AoA
    
    
    return initConvergence


def grad_run():

    # Run the AD version of the code, using known solution from fwd run
    result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)
