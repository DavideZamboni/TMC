import numpy as np
import pandas as pd
from json import load
from MDAnalysis.units import constants, convert
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.due import due, Doi
from MDAnalysis.exceptions import NoDataError
from os import listdir, mkdir
import matplotlib.pyplot as plt
import MDAnalysis as m
from statistics import mean, stdev
from pathlib import Path


def add_charges(u, selection, charge_dict, others = {'MET':{'HT1':0.33, 'HT2':0.33, 'HT3':0.33}, 'GLN':{'OT1':-0.67, 'OT2':-0.67}}):
    charges = load(open(charge_dict))
    sele = u.select_atoms(selection)
    uni_charges = []
    sele_charges = []
    for at in u.atoms:
        if at in sele.atoms:
            rn = at.resname
            an = at.name
            found = False
            if rn in others.keys():
                if an in others[rn].keys():
                    found = True
                    sele_charges.append(charge)
                    uni_charges.append(others[rn][an])
            for (at_name, ff_name, charge) in charges[rn]:
                if an == at_name:
                    uni_charges.append(charge)
                    sele_charges.append(charge)
                    found = True
            if not found:
                print(rn, an, charges[rn])
                raise ModuleNotFoundError(f'The atom {an} of residue {rn} has no associated charge in the dictionary {charge_dict}:\ncharges[{rn}] = {charges[rn]}')
        else:
            uni_charges.append(0)
    u.add_TopologyAttr('charge', uni_charges)
    print(f'Charges succesfully added!\nNumber of atoms: {len(u.atoms)}')
    print(f'Total defined charges: {len(uni_charges)}')
    print(f'Number of selected atoms: {len(sele.atoms)}')
    print(f'Total truly defined charges: {len(sele_charges)}')
    return u

def dielett_calc(u, selection, outfile, center = ['geometry', 'mass'], charge_dict = 'default_name_charges.json'):
    # select atoms:
    atomgroup = u.select_atoms(selection)
    results = pd.DataFrame()
    # get number of frames
    nframes = len(u.trajectory)
    print(f'Number of frames to analyze: {nframes}')

    # if charges not defined, try to add or raise error
    if not hasattr(atomgroup, "charges"):
        try:
            u = add_charges(u, selection, charge_dict)
        except:
            raise NoDataError("No charges defined given atomgroup.")
    
    if center in ['geometry', 'mass']:
        if center == 'mass':
            center_data = np.array([atomgroup.center_of_mass() for ts in u.trajectory])
        if center == 'geometry':
            center_data = np.array([atomgroup.center_of_geometry() for ts in u.trajectory])
        print(f'The total dipole of the selected atoms will be computed from their center of ', center) #, ':\n', center_data.shape, center_data)

    data = {'Time (ns)':[], 'Dip_x':[], 'Dip_y':[], 'Dip_z':[],'Magnitude (eA)':[] , f'co{center[0]}_x':[], f'co{center[0]}_y':[], f'co{center[0]}_z':[]}

    if nframes > 1:
        for ts in u.trajectory:
            my_center = center_data[ts.frame]
            positions = atomgroup.positions - my_center
            #print('Positions: ', positions[:5], '\n')
            #print('Charges: ', atomgroup.charges[:5])
            M = np.dot(atomgroup.charges, positions)
            #print('M', M)
            #print('my_center:', my_center)
            data['Time (ns)'].append(ts.time/1000)
            data['Dip_x'].append(M[0])
            data[f'co{center[0]}_x'].append(my_center[0])
            data['Dip_y'].append(M[1])
            data[f'co{center[0]}_y'].append(my_center[1])
            data['Dip_z'].append(M[2])
            data[f'co{center[0]}_z'].append(my_center[2])
            data['Magnitude (eA)'].append(np.linalg.norm(M))
    else:
        my_center = center_data[0]
        positions = atomgroup.positions - my_center
        M = np.dot(atomgroup.charges, positions)
        #print('M', M)
        #print('my_center:', my_center)
        data['Time (ns)'].append(0)
        data['Dip_x'].append(M[0])
        data[f'co{center[0]}_x'].append(my_center[0])
        data['Dip_y'].append(M[1])
        data[f'co{center[0]}_y'].append(my_center[1])
        data['Dip_z'].append(M[2])
        data[f'co{center[0]}_z'].append(my_center[2])
        data['Magnitude (eA)'].append(np.linalg.norm(M))
    
    results = pd.DataFrame.from_dict(data)
    results.to_csv(outfile, header=True, index=False)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    #print(v1, v2)
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    #print(v1_u, v2_u)
    rad_ang = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    #print(rad_ang)
    return rad_ang

def time_smoothing(vals, times, smooth_window = 10):
    vals_std = []
    smt_vals = []
    start = 0
    while start < len(times) - 2:
        myidxs = [i for i in range(len(times)) if i > start and times[start] + smooth_window >= times[i]]
        myvals = [vals[i] for i in myidxs]
        smt_vals.append(mean(myvals))
        vals_std.append(stdev(myvals))
        start += 1
    
    while len(smt_vals) < len(times):
        last = len(smt_vals) - len(times)
        myvals = vals[last:]
        if len(myvals) > 1:
            smt_vals.append(mean(myvals))
            vals_std.append(stdev(myvals))
        else:
            smt_vals.append(mean(myvals))
            vals_std.append(0)
    
    #print(smt_vals[:-5], vals_std[:-5])
    #print(len(smt_vals), len(times), len(vals_std))
    return smt_vals, vals_std

def plot_compare(inpath, selections, variants, outpath, monomers = ['A', 'B'], to_Debyes = False):
    data = {}
    for var in variants:
        for sel in selections.keys():
            data[f'{var} {sel}'] = f'{inpath}/Dipoles_TMC1_{var}_{"_".join(selections[sel].split())}.csv'
    
    #print(data)

    # General plotting params
    small=14
    medium=16
    large=18
    xl=22
    colors = ["black","dodgerblue","red", "forestgreen", "gold", "darkviolet"]
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    plt.figure(figsize=(6.1,5), dpi=600)

    # Angle between monomers'dipole moment comparison
    angles = {}
    times = {}
    for var in vars:
        dipdata = {}
        for mon in monomers:
            mydf = pd.read_csv(data[f'{var} {mon}'], header=0, index_col=None, usecols=['Time (ns)', 'Dip_x', 'Dip_y', 'Dip_z'])
            dipdata[f'{var} {mon}'] = [np.array([x, y, z]) for (t, x, y, z) in mydf.itertuples(index = False)]
            if var not in times.keys():
                times[var] = [i for i in mydf['Time (ns)']]
        for i in range(len(monomers) - 1):
            for j in range(1, len(monomers)):
                if i != j:
                    angles[f'{var} ({monomers[i]}{monomers[j]})'] = []
                    for m in range(len(dipdata[f'{var} {monomers[j]}'])):
                        #print(var, m, i, j, monomers[i], monomers[j])
                        #print(dipdata[f'{var} {monomers[i]}'][m])
                        #print(dipdata[f'{var} {monomers[j]}'][m])
                        deg_ang = np.degrees(angle_between(dipdata[f'{var} {monomers[i]}'][m], dipdata[f'{var} {monomers[j]}'][m].T))
                        angles[f'{var} ({monomers[i]}{monomers[j]})'].append(deg_ang)
    
    fig, ax=plt.subplots()
    i = 0
    for k in angles.keys():
        y, std = time_smoothing([i for i in angles[k]], [i for i in times[k.split()[0]]], 15)
        ax.plot(times[k.split()[0]], y, label=k, color=colors[i])
        plt.fill_between(times[k.split()[0]], [y[k]-std[k] for k in range(len(y))], [y[k]+std[k] for k in range(len(y))], color=colors[i], alpha=0.5)
        i += 1

    plt.ylabel('Angle (A B, °)')
    plt.xlabel('Time (ns)')
    #plt.xlim((0,145))
    #plt.ylim((0,8))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.minorticks_on()
    #plt.title('Angle between dipoles of each monomer',fontsize = large)
    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    legend.get_frame().set_facecolor("none")
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{outpath}dipoles_angle_smoothed.png", dpi=600)
    plt.close()
    
    if to_Debyes:
        ylab = ' (D)'
        fname = f"{outpath}Debyes_"
    else:
        ylab = ' (e' + r'$\AA$' + ')'
        fname = f"{outpath}eA_"

    # Comparison of magnitude of single monomers
    fig, ax=plt.subplots()
    i = 0
    for var in vars:
        for mon in monomers:
            label=f'{var} {mon}'
            mydf = pd.read_csv(data[f'{var} {mon}'], header=0, index_col=None, usecols=['Time (ns)', 'Magnitude (eA)'])
            if to_Debyes:
                mydf['Magnitude (D)'] = mydf['Magnitude (eA)'] * 4.803
                y, std = time_smoothing([i for i in mydf['Magnitude (D)']], [i for i in mydf['Time (ns)']], 10)
            else:
                y, std = time_smoothing([i for i in mydf['Magnitude (eA)']], [i for i in mydf['Time (ns)']], 10)
            ax.plot(mydf['Time (ns)'], y, label=label, color=colors[i])
            plt.fill_between(mydf['Time (ns)'], [y[k]-std[k] for k in range(len(y))], [y[k]+std[k] for k in range(len(y))], color=colors[i], alpha=0.5)
            i += 1

    plt.ylabel('\N{GREEK SMALL LETTER MU} by monomer' + ylab)
    plt.xlabel('Time (ns)')
    #plt.xlim((0,145))
    #plt.ylim((0,8))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.minorticks_on()
    #plt.title('Dipole moments of each monomer',fontsize = large)
    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    legend.get_frame().set_facecolor("none")
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{fname}monomers_dipole_momenta_smoothed.png", dpi=600)
    plt.close()

    fig, ax=plt.subplots()
    i = 0
    for var in vars:
        mydf = pd.read_csv(data[f'{var} '], header=0, index_col=None, usecols=['Time (ns)', 'Magnitude (eA)'])
        if to_Debyes:
            mydf['Magnitude (D)'] = mydf['Magnitude (eA)'] * 4.803
            y, std = time_smoothing([i for i in mydf['Magnitude (D)']], [i for i in mydf['Time (ns)']], 10)
        else:
            y, std = time_smoothing([i for i in mydf['Magnitude (eA)']], [i for i in mydf['Time (ns)']], 10)

        ax.plot(mydf['Time (ns)'], y, label=f'{var}', color=colors[i])
        plt.fill_between(mydf['Time (ns)'], [y[k]-std[k] for k in range(len(y))], [y[k]+std[k] for k in range(len(y))], color=colors[i], alpha=0.5)
        i += 1

    plt.ylabel('\N{GREEK SMALL LETTER MU}' +'(tot),' + ylab)
    plt.xlabel('Time (ns)')
    #plt.xlim((0,145))
    #plt.ylim((0,8))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.minorticks_on()
    #plt.title('Total dipole moment',fontsize = large)
    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    legend.get_frame().set_facecolor("none")
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{fname}dimer_dipole_moment_smoothed.png", dpi=600)
    plt.close()


basepath = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/to_load/'
outfolder = f'{basepath}/dipoles/'
try:
    mkdir(outfolder)
except OSError:
    True
    
vars = ['wt', 'M654V'] #  
monomers = ['A', 'B']
charge_dict = f'{basepath}scripts/default_name_charges.json'
selections = ['chainID A', 'chainID B', 'protein']
center = 'mass'

for var in vars:
    s = f'{basepath}/xtcs/{var}/{var}_p_l.pdb'
    t = f'{basepath}/xtcs/{var}/{var}_p_l.xtc'
    for sele in selections:
        #u = m.Universe(s, t)
        outfile = f'{outfolder}/Dipoles_TMC1_{var}_{"_".join(sele.split())}.csv'
        #u = add_charges(u, sele, charge_dict, others = {'MET':{'HT1':0.33, 'HT2':0.33, 'HT3':0.33}, 'GLN':{'OT1':-0.67, 'OT2':-0.67}})
        #dielett_calc(u, sele, outfile, center = center)

selections = {'A':'chainID A', 'B':'chainID B', '':'protein'}

plot_compare(outfolder, selections, vars, outfolder, to_Debyes=True)

