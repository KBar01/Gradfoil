import numpy as np
import subprocess
import os

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