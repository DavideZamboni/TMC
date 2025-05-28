#!/usr/bin/env python

import MDAnalysis as mda
import pandas as pd
import os
from os import listdir, mkdir
import matplotlib.pyplot as plt
import matplotlib as mpl
from os.path import isfile, join
import numpy as np
from MDAnalysis.analysis import hole2
import seaborn as sns
from json import dump, load
from statistics import mean, stdev

# Function required for HOLE input
def replace_names(file, 
                  rep = {'C210':'C40', 'C211':'C41', 'C212':'C42', 'C213':'C43', 'C214':'C44', 'C215':'C45', 'C216':'C46', 'C217':'C47', 'C218':'C48',
                               'C310':'C50', 'C311':'C51', 'C312':'C52', 'C313':'C53', 'C314':'C54', 'C315':'C55', 'C316':'C56'},
                  res_chain='POPCX'):
    '''
    HOLE DOES NOT READ ATOM NAMES OF CARBONS LONGER THAN 3 CHARACTER
    USE DICTIONARY REP TO INDICATE ATOM NAMES TO SUBSTITUTE
    USE RES_CHAIN TO IDENTIFY WHICH RESIDUE REQUIRES THIS OPERATION
    '''
    tl = []
    with open(file, 'r') as f:
        ol = f.readlines()
        for line in ol:
            if res_chain in line:
                an = line[12:16]
                if an in rep.keys():
                    newline = line.replace(an, f' {rep[an]}')
                    tl.append(newline)
                else:
                    tl.append(line)
            else:
                tl.append(line)
    f.close()

    with open(file, 'w') as f:
        for l in tl:
            f.write(l)
    f.close()

#### functions for the analysis of HOLE results 
## Functions for obtaining results in .csv format
# From each of the HOLE outputs, retrives center's coordinates, closest atoms (name, resname, resid) and radius (A) and writes it to a .csv file
def get_hole_time_csvs(onlyfiles, trjpath, res_dir, csv_dir, mon):
    '''
    1st step in postprocessing.
    From HOLE outputs, retrives center's coordinates, closest atoms (name, resname, resid) and radius (A).
    Creates for each output file a csv file containing these infos.
    '''
    print('Collecting HOLE results in csvs!')
    #print('usin files: ', onlyfiles)
    ftdic = {}
    df = pd.read_csv(f'{trjpath}/xtcs/Time_frames.csv', header=0, index_col=None)
    for (frame, time) in df.itertuples(index=False):
        ftdic[frame] = time

        
    for f in onlyfiles:
        #print(f)
        frame = int(f[9:-9])
        time = ftdic[frame]
        my_centers = {'Center x':[], 'Center y':[], 'Center z':[], r'Radius $\AA$':[], 'Atom Name':[], 'Resname':[], 'Resid':[]}
        with open(res_dir+f, 'r') as fid:
            ll = fid.readlines()
            #print(ll)
            for l in ll:
                if l.startswith('  at point'):
                    s = l.split()
                    cen_x = float(s[2])
                    cen_y = float(s[3])
                    cen_z = float(s[4])
                    my_centers['Center x'].append(cen_x)
                    my_centers['Center y'].append(cen_y)
                    my_centers['Center z'].append(cen_z)
                elif l.startswith('  closest atom surface'):
                    s = l.split()
                    rad = float(s[3])
                    name = str(s[4])
                    if len(name) > 3:
                        name = str(s[4][:-3])
                        resname = str(s[4][-3:])
                        resid =  int(s[-1])
                    else:
                        resname = s[5]
                        #if len(resname) > 3:
                        #    print('resname' ,resname)
                        resid = int(s[-1])

                    my_centers[r'Radius $\AA$'].append(rad)
                    my_centers['Atom Name'].append(name)
                    my_centers['Resname'].append(resname)
                    my_centers['Resid'].append(resid)
                elif l.startswith('  2nd closest surface'):
                    s = l.split()
                    rad = float(s[3])
                    name = s[4]
                    if len(name)>3:
                        resname = name[-3:]
                        name = name[:4]
                    else:
                        resname = s[5]
                    #print(name, resname)
                    resid =  s[-1]
                    my_centers['Center x'].append(cen_x)
                    my_centers['Center y'].append(cen_y)
                    my_centers['Center z'].append(cen_z)
                    my_centers[r'Radius $\AA$'].append(rad)
                    my_centers['Atom Name'].append(name)
                    my_centers['Resname'].append(resname)
                    my_centers['Resid'].append(resid)
        df = pd.DataFrame.from_dict(my_centers)
        monomers = ['A', 'B']

        df.to_csv(csv_dir + f'monomer_{monomers[mon-1]}_time_{time}ns.csv')

def resids_freq_csv(csv_dir, monomers, charge_dict, polar_c_threshold = 0.3, charged_threshold = 0.7, times=[], add = {'POPC':{'C40':'C210', 'C41':'C211', 'C42':'C212', 'C43':'C213', 'C44':'C214',
                                                                                                                       'C45':'C215', 'C46':'C216', 'C47':'C217', 'C48':'C218', 
                                                                                                                       'C50':'C310', 'C51':'C311', 'C52':'C312', 'C53':'C313', 'C54':'C314', 'C55':'C315', 'C56':'C316'}}
                    ):                                                                                       

    '''
    Writes a csv with Center z and type of closest atoms (Protein or POPC; Polar, Hydrophobic, Acidic or Basic)
    Generally the charge determines the type of atom:
     - if abs(charge) <= polar_c_threshold the atom is Hydrophobic;
     - if polar_c_threshold <= abs(charge) < charged_threshold the atom is Polar;
     - if abs(charge) >= charged_threshold and chrge > 0 the atom is Basic;
     - if abs(charge) >= charged_threshold and chrge < 0 the atom is Acidic;

    There are few exceptions: H atoms of the amino group at the Nter and in the residues ARG, LYS are considered as basic, O atoms of the carboxyl group at the Cter are considered as acidic and also the C of carboxyl groups in GLU and ASP
    If times is specified, performs this operation only on these frames 
    '''

    add = {'POPC':{'C40':'C210', 'C41':'C211', 'C42':'C212', 'C43':'C213', 'C44':'C214', 'C45':'C215', 'C46':'C216', 'C47':'C217', 'C48':'C218', 'C50':'C310', 'C51':'C311', 'C52':'C312', 'C53':'C313', 'C54':'C314', 'C55':'C315', 'C56':'C316'}}                              
    special_cases = {'MET':['HT1', 'HT2', 'HT3', 0.33], 'GLN':['OT1', 'OT2', -0.67]}
    popc_sub = {'Head':['C12', 'C13', 'C14', 'C15', 'H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C', 'H15A', 'H15B', 'H15C', 'H11A', 'H11B', 'H11C', 'C11', 'O11', 'O12', 'O13', 'O14', 'P', 'C1', 'HA', 'HB', 'C2', 'HS', 'C3', 'HX', 'HY', 'O31', 'C31', 'O32', 'O21', 'C21', 'O22']}


    charges = load(open(charge_dict)) # load charge dict
    for resname in add.keys():
        for newname in add[resname]:
            olddata = [i for i in charges[resname] if i[0] == add[resname][newname]][0]
            newdata = [newname, olddata[1], olddata[2]]

            charges[resname].append(newdata)

    
    calls = set()
    for monomer in monomers:
        resid_types = {'z-coordinate (' + r'$\AA$)':[], 'Protein':[], 'Heads':[], 'Tails':[], 'Polar':[], 'Hydrophobic':[], 'Acidic':[], 'Basic':[]}

        data = pd.DataFrame.from_dict({'z-coordinate (' + r'$\AA$)':[], 'Atom Name':[], 'Resid':[], 'Resname':[]})
        features = pd.read_csv(f'{csv_dir}features_{monomer}_alltimes.csv', header=0, index_col=False, usecols=['z-coordinate (' + r'$\AA$)', 'Time (ns)', 'Atom Name', 'Resid', 'Resname'])
        print('features', features.head())
        myts = [float(i) for i in times]
        for t, tdf in features.groupby('Time (ns)'):
            # When looking at frequencies of atoms as pore-lining along a trajectory, the possibility that one atom is present
            # more than once in HOLE results seems wrong => remove duplicates!

            if len(times) == 0:
                tdf = tdf.drop_duplicates(subset = ['Atom Name', 'Resid', 'Resname']) 
                data = pd.concat([data, tdf[['z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resid', 'Resname']]])
            elif t in myts:
                tdf = tdf.drop_duplicates(subset = ['Atom Name', 'Resid', 'Resname'])
                data = pd.concat([data, tdf[['z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resid', 'Resname']]])

        res_df = data[['z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resname']] # Resid not important to determine charge
        min_z = res_df['z-coordinate (' + r'$\AA$)'].min()
        max_z = res_df['z-coordinate (' + r'$\AA$)'].max()
        #print(min_z, max_z)
        sects = [i for i in range(int(min_z-1), int(max_z+2), 1)]
        basres = {'ARG':['NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22'], 'LYS':['NZ', 'HZ1', 'HZ2', 'HZ3']}
        baslist = [k for k in basres.keys()]
        acres = ['ASP', 'GLU']
        for i in range(len(sects[:-1])):
            vals = res_df.loc[(res_df['z-coordinate (' + r'$\AA$)'] >= sects[i]) & (res_df['z-coordinate (' + r'$\AA$)'] < sects[i+1])][['Atom Name', 'Resname']]
            resid_types['z-coordinate (' + r'$\AA$)'].append((sects[i]+sects[i+1])/2)
            prot = 0
            heads = 0
            tails = 0

            acid = 0
            basic = 0
            polar = 0
            hydrophobic = 0

            for (atom_name, resname) in vals.itertuples(index=False):
                if resname == 'POP':
                    if atom_name in popc_sub['Head']: 
                        heads += 1
                    else:
                        tails += 1
                    resname = 'POPC'
                else:
                    prot += 1

                if resname in special_cases.keys() and atom_name in special_cases[resname]:
                    if atom_name.startswith('HT'):
                        calls.add((atom_name, resname, 'Basic'))
                        basic += 1
                    elif atom_name.startswith('OT'):
                        acid += 1
                        calls.add((atom_name, resname, 'Acidic'))

                else:
                    try:
                        charge = float([c[-1] for c in charges[resname] if atom_name == c[0]][0])
                    except:
                        print('Charge not found: ', resname, atom_name)
                    
                    if resname in baslist and atom_name in basres[resname]:
                        calls.add((atom_name, resname, 'Basic'))
                        basic += 1

                    elif abs(charge) >= polar_c_threshold:
                        if abs(charge) >= charged_threshold:
                            if charge < 0 and resname not in basres:
                                acid += 1
                                calls.add((atom_name, resname, 'Acidic'))
                            elif charge > 0 and resname not in acres:
                                basic += 1
                                calls.add((atom_name, resname, 'Basic'))
                            elif atom_name in ['CG', 'CD'] and resname in acres:
                                acid += 1
                                calls.add((atom_name, resname, 'Acidic'))
                            else:
                                if resname in basres:
                                    calls.add((atom_name, resname, 'Basic'))
                                    basic += 1
                                elif resname in acres:
                                    calls.add((atom_name, resname, 'Acidic'))
                                    acid += 1
                        else:
                            calls.add((atom_name, resname, 'Polar'))
                            polar += 1
                    else:
                        hydrophobic += 1
                        calls.add((atom_name, resname, 'Hydrophobic'))

            resid_types['Protein'].append(prot)
            resid_types['Heads'].append(heads)
            resid_types['Tails'].append(tails)

            resid_types['Polar'].append(polar)
            resid_types['Hydrophobic'].append(hydrophobic)
            resid_types['Acidic'].append(acid)
            resid_types['Basic'].append(basic)

        if len(times) == 0:
            res_df_path=csv_dir + f'closest_intertypes_mon_{monomer}_alltimes.csv'
        else:
            res_df_path=csv_dir + f'closest_intertypes_mon_{monomer}_{times[0]}_to_{times[-1]}.csv'
        
        res_df = pd.DataFrame.from_dict(resid_types)
        res_df.to_csv(res_df_path, header=True, index = False)
    
    calls_dic = {'Atom Name':[], 'Resname':[], 'Call':[]}
    for (atname, resname, call) in calls:
        calls_dic['Atom Name'].append(atname)
        calls_dic['Resname'].append(resname)
        calls_dic['Call'].append(call)
            
    calls_path=csv_dir + f'spotted_intertypes.csv'
    cdf = pd.DataFrame.from_dict(calls_dic)
    cdf.to_csv(calls_path, header=True, index = False)

def get_features_csv(csv_dir, features, monomers, times = []):
    '''
    Writes a csv with all the information for each frame in a single csv
    If times is specified, performs this operation only on these frames 
    '''

    for monomer in monomers:
        if len(times) == 0:
            csvs=[f for f in listdir(csv_dir) if isfile(join(csv_dir, f)) and f.endswith('.csv') and f.startswith(f'monomer_{monomer}') ]
            csvs.sort(key = lambda x: float(x.split('_')[-1].rstrip('ns.csv')))
        else:
            csvs=[]
            for i in times:
                csvs.append(f'monomer_{monomer}_time_{i}ns.csv')
        c = 0
        df = pd.DataFrame(columns = features)
        if len(times) != 0:
            for ts in times:
                try:
                    tmp_df = pd.read_csv(csv_dir + f'monomer_{monomer}_time_{ts}ns.csv', header = 0, index_col=0, dtype={'Center x':float, 'Center y':float, 'Center z':float, 'Radius '+ r'$\AA$':float, 'Resid':int, 'Resname':str, 'Atom Name':str})
                    tmp_df['Time (ns)'] = [float(ts) for i in range(len(tmp_df.index))]
                    tmp_df['Monomer'] = [monomer for i in range(len(tmp_df.index))]
                    df = pd.concat([tmp_df, df], axis=0, join='outer', ignore_index=True)
                except FileNotFoundError as error:
                    c += 1
                    print(csv_dir + f'monomer_{monomer}_time_{ts}ns.csv NOT FOUND')
            cen_df_path=csv_dir + f'features_{monomer}_{times[0]}_to_{times[-1]}ns.csv'
            print(f'Total number of files not found: {c}')
        else:
            for x in csvs:
                time =  float(x.split('_')[-1].rstrip('ns.csv'))
                tmp_df = pd.read_csv(csv_dir + x, header = 0, index_col=0, dtype={'Center x':float, 'Center y':float, 'Center z':float, 'Radius '+ r'$\AA$':float, 'Resid':int, 'Resname':str, 'Atom Name':str})
                tmp_df['Time (ns)'] = [time for i in range(len(tmp_df.index))]
                tmp_df['Monomer'] = [monomer for i in range(len(tmp_df.index))]
                df = pd.concat([tmp_df, df], axis=0, join='outer', ignore_index=True)

            cen_df_path=csv_dir + f'features_{monomer}_alltimes.csv'
        df.rename(columns = {'Center x':'x-coordinate (' + r'$\AA$)', 'Center y':'y-coordinate (' + r'$\AA$)', 'Center z': 'z-coordinate (' + r'$\AA$)'}, inplace = True)
        df.to_csv(cen_df_path, header=True, index = False)

def get_obs_csvs(csv_dir, charge_dict, monomers = ['A', 'B'], times=[], obs = ['centers', 'radii', 'resids', 'minz', 'freq', 'features']):
    '''
    Method that retrives the informations required for plotting
    '''

    if obs == 'All_features':
        obs = ['centers', 'features', 'radii', 'resids', 'minz', 'freq']

    print('Obtaining csvs of the observables')
    
    if 'features' in obs:
        print('Doing features')
        features = ['Time (ns)', 'Monomer', 'Center x', 'Center y', 'Center z', 'Radius '+ r'$\AA$', 'Resid', 'Resname', 'Atom Name']
        #get_features_csv(csv_dir, features, monomers, times)
    if 'freq' in obs:
        print('Doing freq')
        resids_freq_csv(csv_dir, monomers, charge_dict, times = times)
    return obs

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

## Plotting Functions
# center densities
def cen_2D_densities(csv_dir, plot_dir, times, monomers = ['A', 'B'], area_limits = (82.5, 117.5), labels = {'z' :[['Cytosol', 0.01, 0.1], ['Exterior', 0.01, 0.9]]}, xlim = (0, 150), ylim = (0, 150), zlim = (0, 200), around_membrane = (75.0, 125.0), step = 10, center = True):
    '''
    Plots density maps with marginal distributions of the points identified as pore centers for both
    the monomers projected on the x-y, x-z and y-z planes
    '''
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=20
    large=25
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    reduced = False

    if len(times) == 0:
        suf = '_alltimes.csv'
    else:
        suf = f'_{times[0]}_to_{times[-1]}.csv'
    cen_df = pd.DataFrame.from_dict({'x-coordinate (' + r'$\AA$)':[], 'y-coordinate (' + r'$\AA$)':[], 'z-coordinate (' + r'$\AA$)':[], 'Monomer':[]})
    for mon in monomers:
        tmp = pd.read_csv(csv_dir + f'features_{mon}{suf}', header=0, index_col=False, usecols=['x-coordinate (' + r'$\AA$)', 'y-coordinate (' + r'$\AA$)', 'z-coordinate (' + r'$\AA$)', 'Monomer'])
        cen_df = pd.concat([cen_df, tmp], ignore_index = True)

    #print(cen_df.shape)
    cen_df = cen_df[(cen_df['x-coordinate (' + r'$\AA$)'] >= xlim[0]) & (cen_df['x-coordinate (' + r'$\AA$)'] <= xlim[1])]
    #print(cen_df.shape)
    cen_df = cen_df[(cen_df['y-coordinate (' + r'$\AA$)'] >= ylim[0]) & (cen_df['y-coordinate (' + r'$\AA$)'] <= ylim[1])]
    #print(cen_df.shape)
    cen_df = cen_df[(cen_df['z-coordinate (' + r'$\AA$)'] >= zlim[0]) & (cen_df['z-coordinate (' + r'$\AA$)'] <= zlim[1])]
    #print(cen_df.shape)

    if len(around_membrane) == 2:
        zlim = (around_membrane[0], around_membrane[1])
        #print(zlim)
        cen_df = cen_df[(cen_df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (cen_df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
        reduced = True

    if center:
        #print(cen_df.describe(), cen_df.head())

        if zlim != None:
            delta = ((zlim[1] - zlim[0]) / 2) + zlim[0]
            zlim = (zlim[0] - delta, zlim[1] - delta)
            #print(zlim)
        else:
            delta = (cen_df['z-coordinate (' + r'$\AA$)'].max() - cen_df['z-coordinate (' + r'$\AA$)'].min()) / 2
            zlim = (cen_df['z-coordinate (' + r'$\AA$)'].min() - delta, cen_df['z-coordinate (' + r'$\AA$)'].max() - delta)
        
        area_limits = (area_limits[0] - delta, area_limits[1] - delta)
        cen_df['z-coordinate (' + r'$\AA$)'] = cen_df['z-coordinate (' + r'$\AA$)'] - delta
        #print(zlim)
        
    if len(times) == 0:
        pname = plot_dir + '2D_centers_density_alltimes'
    else:
        pname = plot_dir + f'2D_centers_density_{times[0]}_to_{times[-1]}'
    
    if reduced:
        pname += f'_z{around_membrane[0]}_to_z{around_membrane[1]}'
    #print('2D cendens xy xlim: ', xlim)
    #print(np.min(cen_df['x-coordinate (' + r'$\AA$)']))
    #print(np.max(cen_df['x-coordinate (' + r'$\AA$)']))
    plot = sns.jointplot(data =cen_df, x=cen_df['x-coordinate (' + r'$\AA$)'], y=cen_df['y-coordinate (' + r'$\AA$)'],  kind='kde', hue='Monomer', thresh = 0.1, fill=True, bw_adjust=.5)
    leg = plt.legend()
    leg.get_frame().set_linewidth(0.0)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plot.ax_joint.set_xticks([round(i) for i in np.linspace(xlim[0], xlim[1], 8)]) # era tutto np.arange(ylim[0], ylim[1] + 1, 25)
    plot.ax_joint.set_xticklabels([round(i) for i in np.linspace(xlim[0], xlim[1], 8)])
    plot.ax_joint.set_yticks([round(i) for i in np.linspace(ylim[0], ylim[1], 8)])
    plot.ax_joint.set_yticklabels([round(i) for i in np.linspace(ylim[0], ylim[1], 8)])

    #plt.tight_layout()
    if labels != None:
        for k in labels.keys():
            if k in ['y', 'x']:
                ax = plot.ax_joint
                for [lab, x_norm, y_norm] in labels[k]:
                    label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (zlim[1] - zlim[0]) + zlim[0]]
                    # print(label)
                    ax.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(pname + '_XY.png', dpi = 600)
    plt.close()
    
    plot = sns.jointplot(data =cen_df, x=cen_df['x-coordinate (' + r'$\AA$)'], y=cen_df['z-coordinate (' + r'$\AA$)'],  kind='kde', hue='Monomer', thresh = 0.1, fill=True, bw_adjust=.5)
    
    leg = plt.legend()
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(xlim)
    plt.ylim(zlim)
    plot.ax_joint.axhline(y = area_limits[0], color='gray', linestyle='--', zorder = 0)
    plot.ax_joint.axhline(y = area_limits[1], color='gray', linestyle='--', zorder = 0)
    plot.ax_joint.set_xticks([round(i) for i in np.linspace(xlim[0], xlim[1], 8)])
    plot.ax_joint.set_xticklabels([round(i) for i in np.linspace(xlim[0], xlim[1], 8)])

    plot.ax_joint.set_yticks([round(i) for i in np.arange(zlim[0], zlim[1] + 1, step)])
    plot.ax_joint.set_yticklabels([round(i) for i in np.arange(zlim[0], zlim[1] + 1, step)])
    plt.tight_layout()
    if labels != None:
        for k in labels.keys():
            if k in ['x', 'z']:
                ax = plot.ax_joint
                for [lab, x_norm, y_norm] in labels[k]:
                    #print(lab, x_norm, y_norm)
                    label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (zlim[1] - zlim[0]) + zlim[0]]
                    #print(label)
                    ax.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(pname + '_XZ.png', dpi = 600)
    plt.close()
    
    plot = sns.jointplot(data =cen_df, x=cen_df['y-coordinate (' + r'$\AA$)'], y=cen_df['z-coordinate (' + r'$\AA$)'],  kind='kde', hue='Monomer', thresh = 0.1, fill=True, bw_adjust=.5)
    plt.xlim(ylim)
    plt.ylim(zlim)
    plot.ax_joint.axhline(y = area_limits[0], color='gray', linestyle='--', zorder = 0)
    plot.ax_joint.axhline(y = area_limits[1], color='gray', linestyle='--', zorder = 0)
    plot.ax_joint.set_xticks([round(i) for i in np.linspace(ylim[0], ylim[1], 8)])
    plot.ax_joint.set_xticklabels([round(i) for i in np.linspace(ylim[0], ylim[1], 8)])
    plot.ax_joint.set_yticks([round(i) for i in np.arange(zlim[0], zlim[1] + 1, step)])
    plot.ax_joint.set_yticklabels([round(i) for i in np.arange(zlim[0], zlim[1] + 1, step)])
    leg = plt.legend()
    leg.get_frame().set_linewidth(0.0)  
    plt.tight_layout()
    if labels != None:
        for k in labels.keys():
            if k in ['y', 'z']:
                ax = plot.ax_joint
                for [lab, x_norm, y_norm] in labels[k]:
                    label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (zlim[1] - zlim[0]) + zlim[0]]
                    #print(label)
                    ax.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(pname + '_YZ.png', dpi = 600)
    plt.close()

def minrad_vs_t(csv_dir, plot_dir, monomers = ['A', 'B'], times = [], around_membrane = (75.0, 125.0), smooth_window = 10):
    '''
    Does line plot of the minimum radius vs t
    '''
    #row = df[df['Radius '+ r'$\AA$'] == df['Radius '+ r'$\AA$'].min()]
    reduced = False
    rad_df = pd.DataFrame.from_dict({'Time (ns)':[], 'z-coordinate (' + r'$\AA$)':[],'Radius '+ r'$\AA$':[], 'Atom Name':[], 'Resid':[], 'Resname':[], 'Monomer':[]})
    if len(times) == 0:
        suf = '_alltimes.csv'
    else:
        suf = f'_{times[0]}_to_{times[-1]}.csv'
    for mon in monomers:
        tmp = pd.read_csv(csv_dir + f'features_{mon}{suf}', header=0, index_col=False, usecols=['z-coordinate (' + r'$\AA$)', 'Monomer', 'Radius '+ r'$\AA$', 'Atom Name', 'Resid', 'Resname', 'Time (ns)'])
        #print(tmp.shape)

        if len(around_membrane) == 2:
            tmp = tmp[(tmp['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (tmp['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
            reduced = True
        #print(tmp.shape)

        for t, mp in tmp.groupby(by='Time (ns)'):
            row = mp[mp['Radius '+ r'$\AA$'] == mp['Radius '+ r'$\AA$'].min()]
            row['Time (ns)'] = [t for i in range(row.shape[0])]
            rad_df = pd.concat([rad_df, row], ignore_index = True)

    #print(rad_df.head())
    gdf = rad_df.groupby('Monomer')
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=16
    large=18
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    i = 0
    cols = [c for c in rad_df.columns]
    for mon, df in gdf:
        df = df.sort_values('Time (ns)')
        #print(df.head(), df.shape)
        if 'Smooth Radius '+ r'$\AA$' not in cols:
            df['Smooth Radius '+ r'$\AA$'], df['Radius std'] = time_smoothing(df['Radius '+ r'$\AA$'].to_list(), df['Time (ns)'].to_list(), smooth_window)
        plt.plot(df['Time (ns)'], df['Smooth Radius '+ r'$\AA$'], label = mon, color=colors[i])
        plt.fill_between(df['Time (ns)'], df['Smooth Radius '+ r'$\AA$']-df['Radius std'], df['Smooth Radius '+ r'$\AA$']+df['Radius std'], color=colors[i], alpha=0.5)
        i += 1
    plt.xlabel('Time (ns)')
    plt.ylabel('Radius '+ r'$\AA$')
    leg = plt.legend()
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(0, 1000)
    plt.ylim(0, 2)
    plt.xticks(range(0, 1001, 200))
    plt.yticks(np.linspace(0, 2, 9))
    plt.tight_layout()

    if not reduced:
        plt.savefig(plot_dir + 'minrad_vs_time.png', dpi=600)
    else:
        plt.savefig(plot_dir + f'minrad_vs_time_z{around_membrane[0]}_to_z{around_membrane[1]}.png', dpi=600)

    plt.close()

def get_distribs(data, title, plot_dir, hue = 'hue', truesmall = None, legend_order = None, xlim = (0, 200), xticks= [i for i in range(0, 201, 25)],
                  area_limits = (82.5, 117.5), labels = [['Cytosol', 0.01, 0.9], ['Exterior', 0.85, 0.9]], around = True, density = True, hist = True, center = True):
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=16
    large=18
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)

    if truesmall is not None:
        plt.rc('legend', fontsize=truesmall)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False

    if center:
        delta = ((xlim[1] - xlim[0]) / 2) + xlim[0]
        #print(xlim, delta)
        xlim = (xlim[0] - delta, xlim[1] - delta)
        #print(xlim)
        area_limits = (area_limits[0] - delta, area_limits[1] - delta)
        # print(xlim)        area_limits = (area_limits[0] - delta, area_limits[1] - delta)
        data['z-coordinate (' + r'$\AA$)'] = data['z-coordinate (' + r'$\AA$)'] - delta
        xticks = [i - delta for i in xticks]        

    if around:
        data = data[(data['z-coordinate (' + r'$\AA$)'] >= xlim[0]) & (data['z-coordinate (' + r'$\AA$)'] <= xlim[1])]

    if density:
        print(title, title.split())
        if title.split()[-4] == 'alltimes':
            bwa = .5
        else:
            bwa = .05

        #if type(times)
        if legend_order is None:
            g = sns.displot(data, x = 'z-coordinate (' + r'$\AA$)', hue = hue, kind="kde", bw_adjust=bwa)
        else:
            g = sns.displot(data, x = 'z-coordinate (' + r'$\AA$)', hue = hue, kind="kde", bw_adjust=bwa, legend=True, hue_order=legend_order)

        plt.ylabel('Density')
        plt.xlabel('z-coordinate (' + r'$\AA$' + ')')
        g.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1)
        plt.xlim(xlim)
        plt.xticks(xticks)
        ys =  plt.ylim()
        # print(ys)
        ylim = (ys[0], ys[1])
        plt.yticks([])
        plt.axvspan(area_limits[0], area_limits[1], color='gray', alpha = 0.5, zorder = 0)
        #plt.tight_layout()
    
    elif hist:
        if legend_order is None:
            g = sns.displot(data, x = 'z-coordinate (' + r'$\AA$)', hue = hue, kind="hist", bins=50)
        else:
            g = sns.displot(data, x = 'z-coordinate (' + r'$\AA$)', hue = hue, kind="hist", bins=50, legend=True, hue_order=legend_order)

        plt.ylabel('Frequency')
        plt.xlabel('z-coordinate (' + r'$\AA$' + ')')
        g.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1)
        plt.xlim(xlim)
        plt.xticks(xticks)
        ys =  plt.ylim()
        # print(ys)
        ylim = (ys[0], ys[1])
        plt.yticks([])
        plt.axvspan(area_limits[0], area_limits[1], color='gray', alpha = 0.5, zorder = 0)

    if labels != None:
        for [lab, x_norm, y_norm] in labels:
            label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (ylim[1] - ylim[0]) + ylim[0]]
            # print(label)
            plt.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(plot_dir + title.replace(' ', '_')+ '.png', dpi=600)
    plt.close()

def resids_freq_plots(csv_dir, monomers, plot_dir, charge_dict, polar_c_threshold = 0.3, charged_threshold = 0.7, xlim = (0, 200), area_limits = (82.5, 117.5), labels = [['Cytosol', 0.01, 0.9], ['Exterior', 0.85, 0.9]], 
                      around_membrane = (75.0, 125.0), times=[], step = 10, add = {'POPC':{'C40':'C210', 'C41':'C211', 'C42':'C212', 'C43':'C213', 'C44':'C214',
                                                                                                                       'C45':'C215', 'C46':'C216', 'C47':'C217', 'C48':'C218', 
                                                                                                                       'C50':'C310', 'C51':'C311', 'C52':'C312', 'C53':'C313', 'C54':'C314', 'C55':'C315', 'C56':'C316'}}, center = True):
    '''
    Plots histogram of frequencies of the types of atoms closest to the pore center vs z-axis coordinate, 2 histograms are produced:
    1 - Protein or POPC vs z-axis
    2 - Acidic, Basic, Polar or Hydrophobic vs z-axis 
    FAI CHARGE 
    '''
    all_single = False
    reduced= False
    aver_comp = {'z-coordinate (' + r'$\AA$)':[], 'hue':[]}
    aver_electype = {'z-coordinate (' + r'$\AA$)':[], 'hue':[]}
    for monomer in monomers:
        if len(times) == 0:
            res_df_path=csv_dir + f'closest_intertypes_mon_{monomer}_alltimes.csv'
            p1 = f'Density of residue closest to the pore center monomer {monomer} alltimes'
            p2 = f'Density of type of residue closest to center monomer {monomer} alltimes'
            p1t = f'Density of residue closest to the pore center whole protein alltimes'
            p2t = f'Density of type of residue closest to center  whole protein alltimes'

        
        elif times == 'all_single':
            all_single = True

        else:
            res_df_path=csv_dir + f'closest_intertypes_mon_{monomer}_{times[0]}_to_{times[-1]}.csv'
            p1 = f'Density of residue closest to the pore center monomer {monomer} {times[0]} to {times[-1]}'
            p2 = f'Density of type of residue closest to center monomer {monomer} {times[0]} to {times[-1]}'
            p1t = f'Density of residue closest to the pore center whole protein {times[0]} to {times[-1]}'
            p2t = f'Density of type of residue closest to center  whole protein {times[0]} to {times[-1]}'


        if len(around_membrane) == 2:
            p1 += f' z{around_membrane[0]} to z{around_membrane[1]}'
            p2 += f' z{around_membrane[0]} to z{around_membrane[1]}'
            p1t += f' z{around_membrane[0]} to z{around_membrane[1]}'
            p2t += f' z{around_membrane[0]} to z{around_membrane[1]}'
            xlim = (around_membrane[0], around_membrane[1])
            xticks = [i for i in np.arange(around_membrane[0], around_membrane[1] + 1, step)]
        else:
            xticks = [i for i in np.arange(xlim[0], xlim[1] + 1, 25)]

        #print(xlim, xticks)
        if all_single:
            #print('HERE!')
            try:
                mkdir(f'{plot_dir}single_frames_freq/')
            except OSError:
                True
            alltimes = [float(f.lstrip(f'monomer_{monomer}_time_').rstrip('ns.csv')) for f in listdir(csv_dir) if f.startswith(f'monomer_{monomer}_time_')] 
            for t in alltimes:
                #print(t)
                rn = f'closest_intertypes_mon_{monomer}_{t}_to_{t}.csv'
                resids_freq_csv(csv_dir, monomers, charge_dict, polar_c_threshold = polar_c_threshold, charged_threshold = charged_threshold, times=[t], add = add)

                res_df_path=csv_dir + rn
                p2 = f'Histogram of type of residue closest to center monomer {monomer} {t} to {t}'
                p2 = f'Density of type of residue closest to center monomer {monomer} {t} to {t}'
                res_df = pd.read_csv(res_df_path, header=0)
                
                if len(around_membrane) == 2:
                    res_df = res_df[(res_df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (res_df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]

                adf = res_df[['z-coordinate (' + r'$\AA$)', 'Acidic']]
                adf['hue'] = ['Acidic' for k in range(adf['Acidic'].shape[0])]
                adf = adf.rename(columns={"Acidic": "Weight"})
                bdf = res_df[['z-coordinate (' + r'$\AA$)', 'Basic']]
                bdf['hue'] = ['Basic' for k in range(bdf['Basic'].shape[0])]
                bdf = bdf.rename(columns={"Basic": "Weight"})
                pdf = res_df[['z-coordinate (' + r'$\AA$)', 'Polar']]
                pdf['hue'] = ['Polar' for k in range(pdf['Polar'].shape[0])]
                pdf = pdf.rename(columns={"Polar": "Weight"})
                hdf = res_df[['z-coordinate (' + r'$\AA$)', 'Hydrophobic']]
                hdf['hue'] = ['Hydrophobic' for k in range(hdf['Hydrophobic'].shape[0])]
                hdf = hdf.rename(columns={"Hydrophobic": "Weight"})
                df = pd.concat([pdf, hdf, adf, bdf])
                print(f'Frequency of type of residue closest to the pore center monomer {monomer} from following dataset:\n', df.shape, df.head())
                data = {'z-coordinate (' + r'$\AA$)':[], 'hue':[]}
                for (z, w, h) in df.itertuples(index=False):
                    for i in range(w):
                        data['z-coordinate (' + r'$\AA$)'].append(z)
                        data['hue'].append(h)
                ddf = pd.DataFrame.from_dict(data)
                get_distribs(ddf, p2, plot_dir+'single_frames_freq/', area_limits = area_limits, xlim = xlim, xticks = xticks, labels = labels, density=False, center = center)
                '''else:
                    print(f'Image for time {t} ns already done!')'''

        else:
            #print(res_df_path)
            res_df = pd.read_csv(res_df_path, header=0)
            #print(res_df.shape)
                            
            if len(around_membrane) == 2:
                res_df = res_df[(res_df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (res_df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
                #print(res_df.shape)
                xlim = around_membrane
            else:
                xlim = (0, 200)

            #print(res_df.head())
            pdf = res_df[['z-coordinate (' + r'$\AA$)', 'Protein']]
            pdf['hue'] = ['Protein' for k in range(pdf['Protein'].shape[0])]
            pdf = pdf.rename(columns={"Protein": "Weight"})

            hdf = res_df[['z-coordinate (' + r'$\AA$)', 'Heads']]
            hdf['hue'] = ['Heads' for k in range(hdf['Heads'].shape[0])]
            hdf = hdf.rename(columns={"Heads": "Weight"})
            pdf = pd.concat([pdf, hdf])

            tdf = res_df[['z-coordinate (' + r'$\AA$)', 'Tails']]
            tdf['hue'] = ['Tails' for k in range(tdf['Tails'].shape[0])]
            tdf = tdf.rename(columns={"Tails": "Weight"})
            pdf = pd.concat([pdf, tdf])

            print(pdf.head())

            data = {'z-coordinate (' + r'$\AA$)':[], 'hue':[]}
            for (z, w, h) in pdf.itertuples(index=False):
                for i in range(w):
                    data['z-coordinate (' + r'$\AA$)'].append(z)
                    data['hue'].append(h)
                    aver_comp['z-coordinate (' + r'$\AA$)'].append(z)
                    aver_comp['hue'].append(h)

            ddf = pd.DataFrame.from_dict(data)
            get_distribs(ddf, p1, plot_dir, xlim = xlim, xticks=xticks, area_limits = area_limits, labels = labels, center = center)

            adf = res_df[['z-coordinate (' + r'$\AA$)', 'Acidic']]
            adf['hue'] = ['Acidic' for k in range(adf['Acidic'].shape[0])]
            adf = adf.rename(columns={"Acidic": "Weight"})

            bdf = res_df[['z-coordinate (' + r'$\AA$)', 'Basic']]
            bdf['hue'] = ['Basic' for k in range(bdf['Basic'].shape[0])]
            bdf = bdf.rename(columns={"Basic": "Weight"})

            pdf = res_df[['z-coordinate (' + r'$\AA$)', 'Polar']]
            pdf['hue'] = ['Polar' for k in range(pdf['Polar'].shape[0])]
            pdf = pdf.rename(columns={"Polar": "Weight"})

            hdf = res_df[['z-coordinate (' + r'$\AA$)', 'Hydrophobic']]
            hdf['hue'] = ['Hydrophobic' for k in range(hdf['Hydrophobic'].shape[0])]
            hdf = hdf.rename(columns={"Hydrophobic": "Weight"})

            df = pd.concat([pdf, hdf, adf, bdf])

            print(f'Frequency of type of residue closest to the pore center monomer {monomer} from following dataset:\n', df.shape, df.head())
            data = {'z-coordinate (' + r'$\AA$)':[], 'hue':[]}
            for (z, w, h) in df.itertuples(index=False):
                for i in range(w):
                    data['z-coordinate (' + r'$\AA$)'].append(z)
                    data['hue'].append(h)
                    aver_electype['z-coordinate (' + r'$\AA$)'].append(z)
                    aver_electype['hue'].append(h)
            ddf = pd.DataFrame.from_dict(data)

            get_distribs(ddf, p2, plot_dir, xlim = xlim, xticks=xticks, area_limits = area_limits, labels = labels, center = center)


    ddf = pd.DataFrame.from_dict(aver_comp)
    get_distribs(ddf, p1t, plot_dir, xlim = xlim, xticks=xticks, area_limits = area_limits, labels = labels, center = center)

    ddf = pd.DataFrame.from_dict(aver_electype)
    get_distribs(ddf, p2t, plot_dir, xlim = xlim, xticks=xticks, area_limits = area_limits, labels = labels, center = center)

def get_charge_profile(csv_dir, monomers, plot_dir, charge_dict, xlim = (0, 200), area_limits = (82.5, 117.5), labels = [['Cytosol', 0.01, 0.9], ['Exterior', 0.85, 0.9]], around_membrane = (70.0, 130.0), times=[], step = 10, center = True):

    add = {'POPC':{'C40':'C210', 'C41':'C211', 'C42':'C212', 'C43':'C213', 'C44':'C214', 'C45':'C215', 'C46':'C216', 'C47':'C217', 'C48':'C218', 'C50':'C310', 'C51':'C311', 'C52':'C312', 'C53':'C313', 'C54':'C314', 'C55':'C315', 'C56':'C316'}}                              
    special_cases = {'MET':['HT1', 'HT2', 'HT3', 0.33], 'GLN':['OT1', 'OT2', -0.67]}
    
    if len(times) == 0:
        fname = 'Average_charge_profile_alltimes.csv'
    else:
        fname = f'Average_charge_profile_{times[0]}_to_{times[-1]}.csv'
    
    res_df_path=csv_dir + fname

    if 'avg_charge_profile.csv' not in listdir(plot_dir):
        charges = load(open(charge_dict)) # load charge dict
        for resname in add.keys():
            for newname in add[resname]:
                olddata = [i for i in charges[resname] if i[0] == add[resname][newname]][0]
                newdata = [newname, olddata[1], olddata[2]]

                charges[resname].append(newdata)

        resid_types = {'z-coordinate (' + r'$\AA$)':[], 'Monomer':[], 'Average charge':[], 'stdev':[], 'N_obs':[]}

        for monomer in monomers:

            data = pd.DataFrame.from_dict({'Time (ns)':[], 'z-coordinate (' + r'$\AA$)':[], 'Atom Name':[], 'Resid':[], 'Resname':[]})
            features = pd.read_csv(f'{csv_dir}features_{monomer}_alltimes.csv', header=0, index_col=False, usecols=['z-coordinate (' + r'$\AA$)', 'Time (ns)', 'Atom Name', 'Resid', 'Resname'])
            for t, tdf in features.groupby('Time (ns)'):
                # When looking at frequencies of atoms as pore-lining along a trajectory, the possibility that one atom is present
                # more than once in HOLE results seems wrong => remove duplicates!

                if len(times) == 0:
                    tdf = tdf.drop_duplicates(subset = ['Atom Name', 'Resid', 'Resname']) 
                    #print(tdf.head())
                    #tdf['Time (ns)'] = [t for i in range(tdf.shape[0])]
                    data = pd.concat([data, tdf[['Time (ns)', 'z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resid', 'Resname']]])

                elif t in times:
                    tdf = tdf.drop_duplicates(subset = ['Atom Name', 'Resid', 'Resname'])
                    #tdf['Time (ns)'] = [t for i in range(tdf.shape[0])]
                    data = pd.concat([data, tdf[['Time (ns)', 'z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resid', 'Resname']]])

            res_df = data[['Time (ns)', 'z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resname']] # Resid not important to determine charge

            min_z = res_df['z-coordinate (' + r'$\AA$)'].min()
            max_z = res_df['z-coordinate (' + r'$\AA$)'].max()

            sects = [i for i in range(int(min_z-1), int(max_z+2), 1)]
            
            for i in range(len(sects[:-1])):
                #print(i, (sects[i] + sects[i + 1])/2)
                vals = res_df.loc[(res_df['z-coordinate (' + r'$\AA$)'] >= sects[i]) & (res_df['z-coordinate (' + r'$\AA$)'] < sects[i+1])][['Time (ns)', 'Atom Name', 'Resname']]
                
                sect_charges = []
                for t, tdf in vals.groupby(by = 'Time (ns)'):
                    time_charges = []
                    for (time, atom_name, resname) in tdf.itertuples(index=False):
                        if resname == 'POP':
                            resname += 'C'
                        if resname in special_cases.keys() and atom_name in special_cases[resname]:
                            time_charges.append(special_cases[resname][-1])
                        else:
                            try:
                                charge = float([c[-1] for c in charges[resname] if atom_name == c[0]][0])
                                time_charges.append(charge)
                            except:
                                print('Charge not found: ', resname, atom_name)

                    if len(time_charges) > 1:
                        sect_charges.append(np.mean(time_charges))
                    elif len(sect_charges) == 1:
                        sect_charges.append(time_charges[0])
                    else:
                        sect_charges.append(0)

                if len(sect_charges) > 1:
                    resid_types['Average charge'].append(np.mean(sect_charges))
                    resid_types['stdev'].append(np.std(sect_charges))
                elif len(sect_charges) == 1:
                    resid_types['Average charge'].append(sect_charges[0])
                    resid_types['stdev'].append(0)
                else:
                    resid_types['Average charge'].append(0)
                    resid_types['stdev'].append(0)

                resid_types['Monomer'].append(monomer)
                resid_types['N_obs'].append(len(sect_charges))
                resid_types['z-coordinate (' + r'$\AA$)'].append((sects[i] + sects[i + 1])/2)

        res_df = pd.DataFrame.from_dict(resid_types)
        res_df.to_csv(res_df_path, header=True, index = False)
        res_df.to_csv(plot_dir + 'avg_charge_profile.csv', header=True, index=False)
    else:
        res_df = pd.read_csv(plot_dir + 'avg_charge_profile.csv', header=0, index_col=False)

    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=16
    large=18
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    i = 0

    if center:
        delta = ((xlim[1] - xlim[0]) / 2) + xlim[0]
        xlim = (xlim[0] - delta, xlim[1] - delta)
        if len(around_membrane) == 2:
            around_membrane = (around_membrane[0] - delta, around_membrane[1] - delta)
            xlim = around_membrane
        
        xticks = np.arange(xlim[0], xlim[1] + 1, step)    
        area_limits = (area_limits[0] - delta, area_limits[1] - delta)

        res_df['z-coordinate (' + r'$\AA$)'] = res_df['z-coordinate (' + r'$\AA$)'] - delta

    #print(xlim, xticks)
    for mon, df in res_df.groupby(by='Monomer'):
        if len(around_membrane) == 2:
            df = df[(df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
            #print(df.shape)
            #print(df.head())

        plt.plot(df['z-coordinate (' + r'$\AA$)'], df['Average charge'], label = mon, color=colors[i])
        plt.fill_between(df['z-coordinate (' + r'$\AA$)'], df['Average charge']-df['stdev'], df['Average charge']+df['stdev'], color=colors[i], alpha=0.5)
        i += 1
    
    plt.xlabel('z-coordinate (' + r'$\AA$)')
    plt.ylabel('Average charge (e)')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(xlim)
    ys =  plt.ylim((-0.4, 0.4))
    # print(ys)
    ylim = (ys[0], ys[1])
    plt.xticks(xticks)
    plt.yticks(np.arange(ys[0], ys[1] + 0.1, 0.2))
    plt.axvspan(area_limits[0], area_limits[1], color='gray', alpha = 0.5, zorder = 0)
    plt.axhline(y = 0, color='black', linestyle='--', zorder = 0)

    plt.tight_layout()

    if labels != None:
        for [lab, x_norm, y_norm] in labels:
            label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (ylim[1] - ylim[0]) + ylim[0]]
            #print(label)
            plt.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(plot_dir + 'charge_profile.png', dpi=600)
    plt.close()
    

def minz_vs_t(csv_dir, plot_dir,  monomers = ['A', 'B'], times = [], area_limits = (82.5, 117.5), xlim = (0, 1000), ylim = (0, 200), labels = [['Cytosol', 0.01, 0.1], ['Exterior', 0.01, 0.9]], around_membrane = (75, 125), step = 10, center = True):
    '''
    Plots:
    1 - z-axis coordinate at which the radius is the minimum for such frame vs time
    2 - Histogram of frequencies of residues closest to the pore center
    '''
    
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=16
    large=18
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False

    minz_df = pd.DataFrame.from_dict({'Time (ns)':[], 'z-coordinate (' + r'$\AA$)':[],'Radius '+ r'$\AA$':[], 'Resid':[], 'Resname':[], 'Monomer':[]})
    if len(times) == 0:
        suf = '_alltimes.csv'
    else:
        suf = f'_{times[0]}_to_{times[-1]}.csv'
    reduced = False
    for mon in monomers:
        tmp = pd.read_csv(csv_dir + f'features_{mon}{suf}', header=0, index_col=False, usecols=['z-coordinate (' + r'$\AA$)', 'Monomer', 'Radius '+ r'$\AA$', 'Resid', 'Resname', 'Time (ns)'])
        if len(around_membrane) == 2:
            tmp = tmp[(tmp['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (tmp['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
            reduced = True

        for t, mp in tmp.groupby(by='Time (ns)'):
            row = mp[mp['Radius '+ r'$\AA$'] == mp['Radius '+ r'$\AA$'].min()]
            row['Time (ns)'] = [t for i in range(row.shape[0])]
            minz_df = pd.concat([minz_df, row], ignore_index = True)

    #print(minz_df.head())

    if center:
        if len(around_membrane) == 2 and ylim != None:
            delta = (ylim[1] - ylim[0]) / 2
            ylim = (around_membrane[0] - delta, around_membrane[1] - delta)
        elif ylim != None:
            delta = (ylim[1] - ylim[0]) / 2
            ylim = (ylim[0] - delta, ylim[1] - delta)
        else:
            delta = (minz_df['z-coordinate (' + r'$\AA$)'].max() - minz_df['z-coordinate (' + r'$\AA$)'].min()) / 2
            ylim = (minz_df['z-coordinate (' + r'$\AA$)'].min() - delta, minz_df['z-coordinate (' + r'$\AA$)'].max() - delta)
        
        area_limits = (area_limits[0] - delta, area_limits[1] - delta)
        minz_df['z-coordinate (' + r'$\AA$)'] = minz_df['z-coordinate (' + r'$\AA$)'] - delta


    gdf = minz_df.groupby('Monomer')
    markers = ['o', '^']
    i = 0
    for mon, df in gdf:
        plt.scatter(df['Time (ns)'], df['z-coordinate (' + r'$\AA$)'], marker = markers[i] ,label = mon, s=1, c=colors[i])
        i += 1
    plt_title = f'Smallest radius z-coord vs time'
    if reduced:
        plt_title += f' z{around_membrane[0]} to z{around_membrane[1]}'
    plt.ylabel('z-coordinate (' + r'$\AA$'+ ')')
    plt.xlabel('Time (ns)')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.axhspan(area_limits[0], area_limits[1], color='gray', alpha = 0.5, zorder = 0)

    plt.axhline(y = area_limits[0], color='gray', linestyle='--', zorder = 0)
    plt.axhline(y = area_limits[1], color='gray', linestyle='--', zorder = 0)
    plt.xticks(np.arange(xlim[0], xlim[1] + 1, 200))
    plt.yticks(np.arange(ylim[0], ylim[1] + 1, step))
    plt.tight_layout()
    
    if labels != None:
        for [lab, x_norm, y_norm] in labels:
            label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (ylim[1] - ylim[0]) + ylim[0]]
            #print(label)
            plt.text(label[1], label[2], label[0], fontsize=small, color='k')

    plt.savefig(plot_dir + plt_title.replace(' ', '_') + '.png', dpi=600)
    plt.close()

    if 'closest residues count'.replace(' ', '_')+ '.csv' not in listdir(plot_dir):
        ress = {}
        minz_df = minz_df.sort_values(by=['Resid'])
        #print(minz_df.head())
        for (t, cen_z, rad, resid, resname, mon) in minz_df.itertuples(index=False):
            if mon not in ress.keys():
                ress[mon] = []
            ress[mon].append(f'{resname} {resid}')
        occurrences = {}
        for mon in ress.keys():
            if mon not in occurrences.keys():
                occurrences[mon] = []
            obs = set(ress[mon])
            for i in obs:
                occurrences[mon].append((i, ress[mon].count(i)))

        freq_data={'Monomer':[], 'Residue':[], 'Count':[]}
        for mon in occurrences.keys():
            myocc = sorted(occurrences[mon], key = lambda x: int(x[1]), reverse=True)
            labels = [i[0] for i in myocc]
            my_ress = [i[1] for i in myocc]
            width = 0.7                 # the width of the bars: can also be len(x) sequence
            plt.bar(labels, my_ress, width)
            plt.ylabel('Frequency')
            plt.xticks(rotation=90, fontsize=5)
            plt_title = f'Frequency of residue indices closest to the pore center monomer {mon}'
            plt.title(plt_title)
            plt.savefig(plot_dir + plt_title.replace(' ', '_')+ '.png', dpi=600)
            plt.close()
            for i in myocc:
                freq_data['Monomer'].append(mon)
                freq_data['Residue'].append(i[0])
                freq_data['Count'].append(i[1])
        df = pd.DataFrame.from_dict(freq_data)
        df.to_csv(plot_dir + 'closest residues count'.replace(' ', '_')+ '.csv', header=True, index=False)

def radius_profile_plotter(csv_dir, plot_dir, monomers, times=[], binsize = 0.5, area_limits = (82.5, 117.5), xlim = (0, 200), ylim = (0, 25), labels = [['Cytosol', 0.01, 0.9], ['Exterior', 0.85, 0.9]], around_membrane = (75, 125), step = 10, center = True):
    '''
    Plots channel profiles averaged over step ns vs z-axis coordinate.
    z-axis binned from minz to maxz in 150 slices
    Channel profiles: mean of radii for such slice + error bar = 1 stdev
    '''
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=14
    medium=16
    large=18
    xl=22
    plt.rc('font',family="arial")
    plt.rc('axes', labelsize=large)
    plt.rc('xtick', labelsize=medium)
    plt.rc('ytick', labelsize=medium)
    plt.rc('legend', fontsize=small)
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
  
    reduced = False
    if center:
        if len(around_membrane) == 2:
            old = [i for i in around_membrane]
            delta = (around_membrane[1] - around_membrane[0]) / 2 + around_membrane[0]
            xlim = (around_membrane[0] - delta, around_membrane[1] - delta)
            ylim = (0, 7)
            yticks = np.arange(ylim[0], ylim[1] + 1, 1)
        else:
            delta = (xlim[1] - xlim[0]) / 2
            xlim = (xlim[0] - delta, xlim[1] - delta)
            yticks = np.arange(0, 25, 2.5)

        area_limits = (area_limits[0] - delta, area_limits[1] - delta)

    xticks = np.arange(xlim[0], xlim[1]+1, step)
    for mon in monomers:
        if len(times) > 0:
            fname = f'pore_profile_{times[0]}_to_{times[-1]}_monomer_{mon}.csv'
        else:
            fname = f'pore_profile_alltimes_monomer_{mon}.csv'

        if fname not in listdir(csv_dir):
            if len(times) > 0:
            
                myfiles = [f for f in listdir(csv_dir) if isfile(join(csv_dir, f)) and f.endswith('ns.csv') and f.startswith(f'monomer_{mon}_time')]
                onlyfiles = []
                for f in myfiles:
                    t = f[15:-6]
                    t = float(t)
                    if t in times:
                        onlyfiles.append(f)
            else:
                onlyfiles = [f for f in listdir(csv_dir) if isfile(join(csv_dir, f)) and f.endswith('ns.csv') and f.startswith(f'monomer_{mon}_time') ]

            #print(f'Pore profile files: ', onlyfiles)
            my_stats=np.ndarray(shape=(0,4))
            for f in onlyfiles:
                mydataset = pd.read_csv(join(csv_dir, f), index_col=0)
                if len(around_membrane) == 2:
                    mydataset = mydataset[(mydataset['Center z'] >= old[0]) & (mydataset['Center z'] <= old[1])]
                #print(mydataset.head())
                x = mydataset.iloc[:,0]
                y = mydataset.iloc[:,1]
                if center:
                    z = mydataset.iloc[:,2] - delta
                else:
                    z = mydataset.iloc[:,2]
                rads = mydataset.iloc[:,3]
                coords = np.array((x, y, z, rads)).T
                my_stats=np.concatenate((my_stats, coords), axis=0)

            data = {'z-coordinate (' + r'$\AA$)':[]}
            # Pore Profile from radii
            profile = my_stats[:, 2:]
            minz = xlim[0]
            maxz = xlim[1]
            #print(minz, maxz)
            mybins = np.arange(minz, maxz + binsize, binsize)
            stats = []
            for i in range(len(mybins)-1):
                data['z-coordinate (' + r'$\AA$)'].append(mybins[i])
                count = []
                for j in profile:
                    if mybins[i]<j[0]<mybins[i+1]:
                        count.append(j[1])
                stats.append((np.mean(count), np.std(count)))
            av_rads = [i[0] for i in stats]
            errors = [i[1] for i in stats]
            data['Radius (' + r'$\AA$)'] = av_rads
            data['std'] = errors
            #print(len(mybins[:-1]), av_rads, len(av_rads), errors, len(errors))
            df = pd.DataFrame.from_dict(data)
            df.to_csv(csv_dir + fname, index=False, header=True)
    
    for mon in monomers:
        if len(times) > 0:
            fname = f'pore_profile_{times[0]}_to_{times[-1]}_monomer_{mon}.csv'
        else:
            fname = f'pore_profile_alltimes_monomer_{mon}.csv'
        mdf = pd.read_csv(csv_dir + fname, header = 0, index_col=False)
        fig, ax = plt.subplots()
        mybins = [i for i in mdf['z-coordinate (' + r'$\AA$)']]
        av_rads = [i for i in mdf['Radius (' + r'$\AA$)']]
        errors = [i for i in mdf['std']]
        ax.bar(mybins, av_rads, yerr=errors, error_kw={'elinewidth':0.3, 'ecolor':'black'})
        plt.xlabel(r'z coordinate ($\AA$)')
        plt.ylabel(r'Pore Radius ($\AA$)')
        #leg = plt.legend()
        #leg.get_frame().set_linewidth(0.0)  
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xticks(xticks)
        plt.yticks(yticks)
        plt.axvspan(area_limits[0], area_limits[1], color='gray', alpha = 0.5, zorder = 0)
        plt.tight_layout()
        if labels != None:
            for [lab, x_norm, y_norm] in labels:
                label = [lab, x_norm * (xlim[1] - xlim[0]) + xlim[0], y_norm * (ylim[1] - ylim[0]) + ylim[0]]
                # print(label)
                plt.text(label[1], label[2], label[0], fontsize=small, color='k')
        if len(times) > 0:
            f'pore_profile_{times[0]}_to_{times[-1]}_{mon}'
        else:
            title = f'pore_profile_whole_trajectory_{mon}'

        if reduced:
            title += f'_z{around_membrane[0]}_to_z{around_membrane[-1]}'

        plt.savefig(plot_dir + title + '.png', dpi=600)
        plt.close()


def get_min_through_time(datapath, trjpath, respath, outpath, var, chargedict, get_csvs = True, plots = True, monomers = ['A', 'B'], obs = 'All_features', times = [], area_limits = (82.5, 117.5), xlim = (0, 140), ylim = (0, 140), zlim = (0, 200), timelim = (0, 1000), smooth_window = 10, labels = [['Cytosol', 0.01, 0.9], ['Exterior', 0.85, 0.9]], around_membrane = (75.0, 125.0), step = 10, center = True):
    print('performing get_min_through_time with times: ', times)
    csv_dir = datapath + f'/csv/{var}/'
    try:
        mkdir(csv_dir)
    except OSError:
        True

    plot_dir = outpath + f'results_{var}/'
    try:
        mkdir(plot_dir)
    except OSError as error:
        True

    if len(around_membrane) == 2:
        plot_dir = outpath + f'results_{var}/reduced/'
        try:
            mkdir(plot_dir)
        except OSError as error:
            True


    if obs == 'All_features':
        obs = ['centers', 'freq', 'radii', 'minz', 'pore_profiles', 'features', 'charge_profile']
    
    if get_csvs:
        for mon in range(1, len(monomers) + 1):
            if len(times) == 0 or type(times) == str:
                onlyfiles = [f for f in listdir(respath) if isfile(join(respath, f)) and f.endswith('_hole.out') and f.startswith(f'Ch{mon}_frame') ]
                #get_hole_time_csvs(onlyfiles, trjpath, respath, csv_dir, mon)
        if type(times) == str:
            get_obs_csvs(csv_dir, chargedict, monomers = ['A', 'B'], times = [], obs = obs)
        else:
            get_obs_csvs(csv_dir, chargedict, monomers = ['A', 'B'], times = times, obs = obs)

    #### POSSIBILITIES 'centers', 'radii', 'minz', 'freq', 'pore_profiles', 'charge_profile'
    
    if plots:
        if 'centers' in obs:
            print('Doing centers')
            cen_2D_densities(csv_dir, plot_dir, times,  area_limits=area_limits, xlim = xlim, ylim = ylim, zlim = zlim, around_membrane = around_membrane, step = step, center = center)

        if 'freq' in obs:
            print('Doing resids freq')
            resids_freq_plots(csv_dir, monomers, plot_dir, chargedict, area_limits=area_limits, labels=labels, around_membrane = around_membrane, times = times, step = step, center = center)

        if 'radii' in obs:
            print('Doing radii')
            minrad_vs_t(csv_dir, plot_dir, monomers=monomers, times=times, around_membrane = around_membrane, smooth_window=smooth_window)

        if 'minz' in obs:
            print('Doing minz')
            minz_vs_t(csv_dir, plot_dir, monomers=monomers,xlim = timelim, ylim=zlim, area_limits=area_limits, around_membrane = around_membrane, step = step, times = times, center = center)

        if 'pore_profiles' in obs:
            print('Doing pore_profiles')
            radius_profile_plotter(csv_dir,  plot_dir, monomers, times = times, xlim = zlim, around_membrane = around_membrane, step = step, labels = labels, center = center)

        if 'charge_profile' in obs:
            print('Doing charge_profile')
            get_charge_profile(csv_dir, monomers, plot_dir, charge_dict, xlim = zlim, area_limits = area_limits, labels = labels, around_membrane = around_membrane, times=times, step = step, center = center)

vars = ['wt', 'M654V']
basepath = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/to_load/'
path = f'{basepath}xtcs/'                   # path where to find trajectories
outpath=f'{basepath}hole/'                      # path where to save all the outfiles
charge_dict = f'{basepath}/scripts/default_name_charges.json' # .json file containing atom charges
channel_helices = ['S4', 'S5', 'S6', 'S7']
monomers = ['A', 'B']
base_selection = 'resname POPC or protein'

## PERFORM ACTUAL HOLE ANALYSIS
# 1 - initialise data acquisition and output
# General HOLE input
cvect = [0.0, 0.0, 1.0] # channel vector, generally aligned to the z-axis
sample = 0.2            # step-size along cvect when performing the analysis (A)
endrad = 22.0           # stop analysis when radius > 22 A

try:
    mkdir(outpath) # General output path
except OSError:
    True

try:
    mkdir(outpath + '/pdbs') # Single-structure pdbs output path, required for HOLE analysis
except OSError as error:
    True
try:
    mkdir(outpath + '/hole_results') # Folder for HOLE results
except OSError as error:
    True

# iterate over variants
for var in vars:
    
    # definition of variant-specific input and output
    inpath= f'{path}{var}/' 
    framesfile = f'{inpath}Newtraj_timframes.csv'
    strpath = outpath + f'/pdbs/{var}/' # Dividing for each of the variants
    try:
        mkdir(strpath)
    except OSError as error:
        True

    respath = outpath + f'/hole_results/{var}/'
    try:
        mkdir(respath)
    except OSError as error:
        True

    # Select input data
    struct = f'{inpath}{var}_p_l.pdb'
    traj = f'{inpath}{var}_p_l.xtc'
    
    # to select the correct atoms of each of the monomers
    monomer_residues = 760
    if var == 'wt':
        index = 12344
    else:
        index = 12343

    ### 2 - HOLE analysis
    os.chdir(respath) # to have outputs in the correct folder
    u = mda.Universe(struct, traj)

    # obtaining database with times and frames, useful for downstream analyses
    if 'Newtraj_timframes.csv' not in os.listdir(inpath):
        newtimerames = {'Time (ns)':[], 'Frame':[]}
        for ts in u.trajectory:
            newtimerames['Frame'].append(ts.frame)
            newtimerames['Time (ns)'].append(ts.time / 1000)

        timeframes_df = pd.DataFrame.from_dict(newtimerames)
        timeframes_df.to_csv(framesfile, header = True, index=False)
        
    timeframes_df = pd.read_csv(framesfile, header = 0, index_col=False)
    times = timeframes_df['Time (ns)'].tolist()
    times.append(0.0)
    times.sort()

    frames = timeframes_df['Frame'].tolist() # perform on all frames
    frames.append(0)
    frames.sort()

    for frame in frames:
        u.trajectory[frame]
        
        if f'TMC1_{var}_frame_{frame}.pdb' not in listdir(strpath): # check if structure already in structure folder, if not create it
            #print('writing structure')
            my_str = u.select_atoms(base_selection)
            my_str.write(f'{strpath}TMC1_{var}_frame_{frame}.pdb')

        # iterate over each of the monomers (1 channel per monomer) and perform HOLE analysis
        c = 0
        for i in range(len(monomers)):
            # Find structure-specific HOLE input data 
            resids = [454, 528, 414, 583]  # Residues in pore-forming helices whose average position will give starting point of HOLE analysis
            if i == 0:
                s = f'protein and index 0:{index} and ('
                chan=1
            elif i > 0:
                s = f'protein and index {index+1}:{index*2} and ('
                chan=2
            for i in resids:
                s += f'resid {i+c} or '
            s = s[:-3] + ')'
            c += monomer_residues

            outfile = f'Ch{chan}_frame{frame}_hole.out'
            #if not isfile(join(f'{respath}', outfile)):
            aas = u.select_atoms(s) # selection of the residues
            # HOLE input                
            cp = aas.center_of_mass() # point within the channel 
            f = f'{strpath}TMC1_{var}_frame_{frame}.pdb'
            # Replace POPC atom names in the .pdb file as HOLE reads atom names longer than 3 characters
            replace_names(f)
            if frame == frames[0]:
                ha = hole2.hole(f, infile=f'{respath}example_IN_Ch_{chan}',  outfile = outfile, executable='/home/davide/hole2/exe/hole', sphpdb_file = f'{respath}Ch{chan}_frame{frame}_hole.sph',cpoint=cp, cvect=cvect, sample=sample, end_radius=endrad, output_level=0, keep_files=True) # 
            else:
                ha = hole2.hole(f, outfile = f'{respath}Ch{chan}_frame{frame}_hole.out', executable='/home/davide/hole2/exe/hole', sphpdb_file = f'{respath}Ch{chan}_frame{frame}_hole.sph',cpoint=cp, cvect=cvect, sample=sample, end_radius=endrad, output_level=0, keep_files=True)
            #else:
            #    print(f'File for frame {frame} already present in {respath}\n{outfile}')
    print(var)

## PERFORM ANALYSIS OF HOLE RESULTS
vars = ['wt','M654V'] #

for var in vars:
    inpath = path + f'/{var}/'
    
    try:
        mkdir(outpath + '/csv/')
    except OSError as error:
        True

    respath = outpath + f'/hole_results/{var}/'
    obs = ['centers', 'radii', 'minz', 'freq', 'pore_profiles', 'charge_profile'] #'centers', 'radii', 'minz', 'freq', 'pore_profiles', 'charge_profile'
    get_min_through_time(outpath, inpath, respath, outpath, var, charge_dict, get_csvs = True, plots = True, monomers = ['A', 'B'], obs = obs, around_membrane = (70, 130), step = 10, center = True)
    
