import json
import subprocess
import numpy as np


coords = np.genfromtxt('n2412coords.txt',delimiter=',')

xcoords = coords[0,:].tolist()
#ycoords = coords[1,:].tolist() 


#ycoords_numpy = np.copy(coords[1,:])
#ycoords = ycoords_numpy.tolist()



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

    # Run the executable for first time, no restarting
    result = subprocess.run(["./build/CFoil_AD"], capture_output=True, text=True)
    #initConvergence = result.returncode


##################### single run  #######################
xcoords = coords[0,:].tolist()
ycoords = coords[1,:].tolist() 
conv = standard_run(3,1e6,0.0,xcoords,ycoords)
grad_run()

############################################ CL alfpha ###########################################



def do_clalf():
    alpha_list = np.linspace(-5,14,int(((14+5)/0.1)+1))

    coords = np.genfromtxt('n2412coords.txt',delimiter=',')

    xcoords = coords[0,:].tolist()
    ycoords = coords[1,:].tolist() 
    n = len(alpha_list)
    print(alpha_list)
    clvals = np.zeros([n,3])
    for i in range(n):

        print("Doing: " +str(alpha_list[i]))
        conv = standard_run(alpha_list[i],1e6,0.0,xcoords,ycoords)


        # Read the JSON file written by C++
        with open("out.json", "r") as f:
            data = json.load(f)

        # Access values
        cl = data["CL"]
        clvals[i,2] = conv 
        clvals[i,1] = cl
        clvals[i,0] = alpha_list[i]

    np.savetxt("Cn2412clalf_05step.csv",clvals,delimiter=',')



#FD Code
def do_FD():
    order = np.logspace(-4,-8,5,base=10)

    print(order)

    gradient = np.zeros([200,len(order)])

    for i in range(len(order)): # loop over step size 
        
        print(order[i])
        for coord in range(200):    #loop over y coordinates
            
            steps = np.zeros(2)
            for FDstep in range(2):
                
                if FDstep==0:
                    nudge = order[i]
                else:
                    nudge = -order[i]
                alpha = 3.0
                Re = 1.0e6
                Ma = 0.0
                
                
                ycoords_numpy = np.copy(coords[1,:])
                ycoords_numpy[coord] += nudge
                ycoords = ycoords_numpy.tolist()
                
                standard_run(alpha,Re,Ma,xcoords,ycoords)

                
                # Read the JSON file written by C++
                with open("out.json", "r") as f:
                    data = json.load(f)

                # Access values
                upperTau = data["upperTauMax"]
                
                steps[FDstep] = upperTau
        
            gradient[coord,i] = (steps[0] - steps[1]) / (2*order[i])

    print(gradient)



    np.savetxt("FDUpperTau.csv",gradient,delimiter=',')


"""
with open("ad_gradients.json", "r") as f:
        data = json.load(f)

AD_taumax = np.array(data["d tauMaxUpper / d ycoords"])

FD_taumax = np.genfromtxt("FDUpperTau.csv",delimiter=',')


for i in range(FD_taumax.shape[1]):
    
    diff = np.max(np.abs((AD_taumax - FD_taumax[:,i])/FD_taumax[:,i]) * 100)
    print(diff)

"""

