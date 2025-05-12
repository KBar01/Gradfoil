import json
import subprocess
import numpy as np

from pkg_resources import resource_filename
pth = resource_filename(__name__, "CFoil_fwd")
coords = np.genfromtxt('n2412coords.txt',delimiter=',')

xcoords = coords[0,:].tolist()
#ycoords = coords[1,:].tolist() 

def standard_run(alphaDeg,Re,Ma,xcoords,ycoords):

    data = {
        "xcoords": xcoords,
        "ycoords": ycoords,
        "alpha_degrees": alphaDeg,
        "Re": Re,
        "Ma": Ma,
        "restart": 0
    }

    # Write JSON input
    with open("input.json", "w") as f:
        json.dump(data, f)

    # Run the executable for first time, no restarting
    result = subprocess.run(["./build/CFoil_fwd"], capture_output=True, text=True)
    print(result)
    initConvergence = result.returncode

    stepsize = 0.5
    if (not initConvergence):
        tempalf = alphaDeg - stepsize ; 
        for i in range(50):
            
            with open("input.json", "r") as f:
                data = json.load(f)

            # Modify only the desired fields
            #data["restart"] = 1  # or 0
            data["alpha_degrees"] = tempalf  # assign your new value here

            # Write back the updated JSON
            with open("input.json", "w") as f:
                json.dump(data, f) 


            result = subprocess.run(["./build/CFoil_fwd"], capture_output=True, text=True)
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

            with open("input.json", "r") as f:
                    data = json.load(f)

            data["restart"] = 1 # using restart fle
            data["alpha_degrees"] = fwdalf  # trying different alfa
            
            with open("input.json", "w") as f:
                json.dump(data, f, indent=4)  # indent is optional, just for readability

            result = subprocess.run(["./build/CFoil_fwd"], capture_output=True, text=True)
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
    result = subprocess.run(["./build/CFoil_AD"], capture_output=True, text=True)
