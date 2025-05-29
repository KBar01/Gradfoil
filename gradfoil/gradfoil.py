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
        "restart": 0,
        "xfoilStart":0,
        "xfoilgetpoints":0
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

###################################################### doing xfoil run #######################################################
def extract_foil_data(array):

    firstSection = []
    i = 0
    xcurr = array[i,0] 
    while xcurr <= 1.0:
        val = array[i,1]
        firstSection.append(val)
        i += 1 
        xcurr = array[i,0]
    
    #skip over wake section
    xprev = array[i-1,0]
    while xcurr > xprev:
        i += 1
        xcurr = array[i,0]
        xprev = array[i-1,0]

    #extract top theta
    secondSection = []
    while xcurr <= 1.0:
        val = array[i,1]
        secondSection.append(val)
        i += 1 
        if i >= array.shape[0]:
            break
        xcurr = array[i,0]

    return np.array(firstSection),np.array(secondSection)


def extract_ncrit_data(array):

    firstSection = []
    i = 0
    nprev = 0
    ncurr = array[i,1]
    while ncurr>=nprev:
        val = array[i,1]
        firstSection.append(val)
        i += 1 
        ncurr = array[i,1]
        nprev = array[i-1,1]


    #extract top theta
    secondSection = []
    while i < array.shape[0]:
        secondSection.append(array[i,1])
        i += 1

    return np.array(firstSection),np.array(secondSection)


def run_xfoil_get_BL_states(coords, alpha, Re, mach):
    """Run XFOIL on a 200-point airfoil coordinate array and extract boundary layer data."""

    # Filenames
    dat_file = 'airfoil.dat'
    xfoil_in = 'xfoil_input.txt'
    files = {
        'n': 'bl_n.txt',
        'ct': 'bl_ct.txt',
        'top': 'bl_top.txt',
        'bot': 'bl_bot.txt',
        'ue': 'bl_ue.txt'
    }

    # Write airfoil.dat file
    with open(dat_file, 'w') as f:
        f.write('Foil\n')
        for x, y in coords:
            f.write(f"{x:.6f} {y:.6f}\n")

    # Write XFOIL input script
    with open(xfoil_in, 'w') as f:
        f.write(f"LOAD {dat_file}\n")
        f.write("OPER\n")
        f.write(f"VISC {Re:.0f}\n")
        f.write(f"MACH {mach:.4f}\n")
        f.write("ITER 100\n")
        f.write(f"ALFA {alpha}\n")
        f.write("VPLO\n")
        f.write("N\nDUMP\n" + files['n'] + "\n")
        f.write("CT\nDUMP\n" + files['ct'] + "\n")
        f.write("DT\nDUMP\n" + files['top'] + "\n")
        f.write("DB\nDUMP\n" + files['bot'] + "\n")
        f.write("UE\nDUMP\n" + files['ue'] + "\n\nQUIT\n")

    # Run XFOIL
    subprocess.run([r'/mnt/c/Users/pa20830/XFOIL6.99/xfoil.exe'],
                   input=open(xfoil_in).read(),
                   text=True)

    def load_column_data(file):
        data = []
        with open(file, 'r') as f:
            for line in f:
                try:
                    a, b = map(float, line.strip().split())
                    data.append((a, b))
                except:
                    continue
        return np.array(data)

    # Load all BL data
    n_data = load_column_data(files['n'])
    ct_data = load_column_data(files['ct'])
    top_data = load_column_data(files['top'])
    bot_data = load_column_data(files['bot'])
    ue_data = load_column_data(files['ue'])

    
    topDelta,topTheta = extract_foil_data(top_data)
    botDelta,botTheta = extract_foil_data(bot_data)

    topN,botN = extract_ncrit_data(n_data)
    topCt,botCt = extract_foil_data(ct_data)
    
    topUe,botUe = extract_foil_data(ue_data)
   
    th = np.concatenate((np.flip(botTheta), topTheta))
    ds = np.concatenate((np.flip(botDelta),topDelta))
    
    sa = np.concatenate((np.flip(botCt), np.flip(botN), topN, topCt))
    ue = np.concatenate((np.flip(botUe), topUe))

    out = np.vstack((th, ds, sa, ue))  # shape (4, N)

    # Create integer mask with 1s for Ct arrays and 0s for N arrays
    turb = np.concatenate((
        np.ones(len(botCt), dtype=int),   # 1s for flipped botCt
        np.zeros(len(botN), dtype=int),   # 0s for flipped botN
        np.zeros(len(topN), dtype=int),   # 0s for topN
        np.ones(len(topCt), dtype=int)    # 1s for topCt
    ))

    # Optionally clean up temporary files
    #for f in [dat_file, xfoil_in] + list(files.values()):
    #    if os.path.exists(f):
    #        os.remove(f)

    return out,turb  # shape (4, N)

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
