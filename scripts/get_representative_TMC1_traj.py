import sys
import os
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align
import numpy as np
import pandas as pd
        
def new_get_representative_frames(my_2d_rmsd_mat, frametimes, step = 10, actstart = 0): #STEP IN #frames, 100 fraems = 1 ns
    '''
    Given a MDAnalysis universe, a dictionary of selection stings (as the ones coming from sel_from_dict), the monomer to analyze (keys of selections)
    the path to working directory and a number representing a window of frames (eg 10 default is to get windows of 1 ns on trajectories with timestep 
    of 10 ps), this script returns, for each consequent window of the trajectory, the number of the frame that has the lowest RMSD with respect to 
    all the other frames.
    '''
    to_analyze={'Time (ns)':[], 'Frame':[]}
    mymat = np.loadtxt(my_2d_rmsd_mat)
    if len(frametimes) == mymat.shape[0]:
        # from 1 as frame 0 will always be starting structure
        if actstart != 0:
            mystart = actstart
        else:
            mystart = 1
        for i in range(mystart, mymat.shape[0], step):
            matvals = mymat[i:i+step, i:i+step]
            rmsd_sum = np.sum(matvals, axis = 0)
            mymin = np.argmin(rmsd_sum)
            to_analyze['Frame'].append(frametimes[i + mymin][0])
            to_analyze['Time (ns)'].append(frametimes[i + mymin][1])

    return to_analyze

def new_get_representative_trj(struc, trj, path, outsel, my_2d_rmsd_mat, idx, outnames, step = 10, monomers = ['A', 'B']):
    outsel = outsel.split('|')
    universe = mda.Universe(struc, trj)
    
    if 'Allframes.csv' not in os.listdir(path):
        frametimes = {'Time (ns)':[], 'Frame':[]}
        for ts in universe.trajectory:
            frametimes['Frame'].append(ts.frame)
            frametimes['Time (ns)'].append(ts.time/1000)
        df = pd.DataFrame.from_dict(frametimes)
        df.to_csv(path + '/Allframes.csv', header=True, index=False)

    
    if 'Repr_frames.csv' not in os.listdir(path):
        df = pd.read_csv(path + '/Allframes.csv', header=0, index_col=False)
        frametimes = df.to_dict(orient='list')
        sorted_frametimes = [(frametimes['Frame'][i], frametimes['Time (ns)'][i]) for i in range(len(frametimes['Frame']))]
        sorted_frametimes.sort(key = lambda x: int(x[0]))

        print(f'Sorted frames and times: {sorted_frametimes}\nEvaluating frames of the trajectory, output in {path}Repr_frames.csv')
        to_analyze = new_get_representative_frames(my_2d_rmsd_mat, sorted_frametimes, step = step)
        df = pd.DataFrame.from_dict(to_analyze)
        df.to_csv(path + '/Repr_frames.csv', header=True, index=False)

    else:
        df = pd.read_csv(path + '/Repr_frames.csv', header=0, index_col=False)
        
    my_f = df['Frame'].to_list()

    sele = {}

    for i in range(len(monomers)):
        sele[monomers[i]] = set([(at.resid, at.resname) for at in universe.select_atoms(f'index {(idx*i)}:{(idx*(i+1) - 1)}').residues])

    chains = []
    for at in universe.atoms:
        Found = False
        for k in sele.keys():
            if (at.resid, at.resname) in sele[k]:
                chains.append(k)
                Found = True
        if not Found:
            chains.append('X')

    universe.add_TopologyAttr('chainID', chains)

    atsels = {}
    for i in outsel:
        out = universe.select_atoms(i)
        atsels[i] = out
    
    my_f.sort()
    for i in range(len(outsel)):
        sel = outsel[i]
        name = outnames[i]
        universe.trajectory[my_f[0]]
        atsels[sel].write(path + f'{name}.gro')
        atsels[sel].write(path + f'{name}.pdb', bonds = 'conect')
        if f'{name}.xtc' not in os.listdir(path):
            atsels[sel].write(path + f'{name}.xtc', frames=universe.trajectory[my_f])
        else:
            atsels[sel].write(path + f'{name}_reduced.xtc', frames=universe.trajectory[my_f])


basepath=str(sys.argv[1])
vars=['wt', 'M654V'] 
for var in vars:
    mons=['A','B']
    basepath='/home/davide/storage/lab_dati/davide_TMC1/definitivi/'
    struc=f'{basepath}trj_data/{var}/TMC1_{var}_npt.tpr'
    trj=f'{basepath}trj_data/{var}/xtcs/TMC1_{var}_noPBC.xtc'
    step=10
    outsel='all|protein or resname POPC or resname POT or resname CLA|name CA'
    outnames=[f'{var}_{i}' for i in  'noPBC p_l CA'.split()]

    path=f'{basepath}/to_load/xtcs/{var}/'
    resnum = 760
    my_2d_rmsd_mat = f'{basepath}trj_data/{var}/xtcs/TMC1_{var}_noPBC_channel_CA.txt'


    print(f'''
    	var: {var}
    	mons: {mons}
    	struc: {struc}
    	trj: {trj}
    	step: {step}
    	outsel: {outsel}
    	outnames: {outnames}
    	path: {path}
          ''')

    if 'wt' in var:
        index = 12344
    elif 'M654V' in var:
        index = 12343

    #new_get_representative_trj(struc, trj, path, outsel, my_2d_rmsd_mat, index, outnames, step = step, monomers = mons)

    framesfile = f'{path}Repr_frames.csv'

    timeframes_df = pd.read_csv(framesfile, header = 0, index_col=False)
    times = timeframes_df['Time (ns)'].tolist()
    times.append(0.0)
    times.sort()
    print(times)
    frames = timeframes_df['Frame'].tolist() # perform on all frames
    frames.append(0)
    frames.sort()
    print(frames)

