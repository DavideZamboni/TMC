#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from os import listdir, mkdir


def densmap_plotter(norm_mat, filename, plane, cm, tick_dict):

    fig, ax = plt.subplots()
    shape = norm_mat.shape
    
    xlabels = tick_dict[plane[0]]
    ylabels = tick_dict[plane[1]]
    xticks = np.linspace(0, shape[1], len(xlabels))
    yticks = np.linspace(0, shape[0], len(ylabels))


    im = ax.imshow(norm_mat, cmap=cm)
    plt.xticks(xticks, xlabels, size=16, rotation=45)
    plt.xlabel(f'{plane[0]}-axis ' + r'$\AA$', size=20)  
    plt.yticks(yticks, ylabels, size=16)
    plt.ylabel(f'{plane[1]}-axis ' + r'$\AA$', size=20)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Density', rotation=-90, va="bottom", fontsize=10)
    plt.tight_layout()
    plt.savefig(filename, dpi=1000)
    plt.close()


def contours_plotter(mats, contours, subject, filename, plane, cm, tick_dict, colors = ['blue', 'green', 'purple', 'red', 'black']):
    # contour of each helix, different colors + black for non pore-forming
    if subject != 'prot_et_mem':
        norm_mat = mats[f'{subject}_{plane}']
    else:
        norm_mat = 1 - mats[f'{subject}_{plane}']

    shape = norm_mat.shape
    
    xlabels = tick_dict[plane[0]]
    ylabels = tick_dict[plane[1]]
    xticks = np.linspace(0, shape[1], len(xlabels))
    yticks = np.linspace(0, shape[0], len(ylabels))
    i = 0

    for c in contours:
        mycontours = plt.contour(mats[f'{c}_{plane}'], levels=[0.2], linewidths=2, colors=colors[i], alpha=0.5)
        plt.clabel(mycontours, inline=True, fontsize=4)
        i += 1

    plt.imshow(norm_mat, cmap=cm, vmin=0)
    plt.xticks(xticks, xlabels, size=16, rotation=45)
    plt.xlabel(f'{plane[0]}-axis ' + r'$\AA$', size=20)
    plt.yticks(yticks, ylabels, size=16)
    plt.ylabel(f'{plane[1]}-axis ' + r'$\AA$', size=20)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(filename, dpi=1000)
    plt.close()

def density_plotter(path, outpath, var, myshape, targets, axes = ['z'], contours = [], densities = False):
    cm = 'cool'
    mats = {}
    poss_axes = ['x', 'y', 'z']
    mykeys = set()
    for i in range(len(axes)):
        axis = axes[i]
        matrices = [path + f'{axis}_axis/' + f for f in listdir(path+f'{axis}_axis/') if f'_{axis}.txt' in f]
        plane = ''.join([k for k in poss_axes if k != axis])
    
        for m in matrices:
            if f'lipids_densmap_{axis}' in m:
                subj = 'lipids'
                mykeys.add('lipids')
            elif f'S4_densmap_{axis}' in m:
                mykeys.add('S4')
                subj = 'S4'
            elif f'S5_densmap_{axis}' in m:
                subj = 'S5'
                mykeys.add('S5')

            elif f'S6_densmap_{axis}' in m:
                subj = 'S6'
                mykeys.add('S6')

            elif f'S7_densmap_{axis}' in m:
                subj = 'S7'
                mykeys.add('S7')

            elif f'no_channel_densmap_{axis}' in m:
                subj = 'non pore-forming TM'
                mykeys.add('non pore-forming TM')

            elif f'POT_densmap_{axis}' in m:
                subj='potassium'
                mykeys.add('potassium')

            elif f'CLA_densmap_{axis}' in m:
                subj='chloride'
                mykeys.add('chloride')

            elif f'protein_densmap_{axis}' in m:
                subj='protein'
                mykeys.add('protein')

            elif f'prot_et_mem_densmap_{axis}' in m:
                subj='prot_et_mem'
                mykeys.add('prot_et_mem')

            density_mat = np.loadtxt(m)
            #print(f'Found density matrix file: {m.split("/")[-1]}, corresponding subject: {subj}\nMatrix shape: {density_mat.shape}')
            mymax = np.max(density_mat)
            mymin = np.min(density_mat)
            #print(f'maxmium value prenorm: {mymax} minimum value prenorm: {mymin}')
            norm_mat = (density_mat - mymin) /(mymax - mymin)
            mymax = np.max(norm_mat)
            mymin = np.min(norm_mat)
            #print(f'maxmium value postnorm: {mymax} minimum value postnorm: {mymin}')
            
            mats[f'{subj}_{plane}'] = norm_mat

    tick_dict = {
                'x':[-70, -50, -30, -10, 10, 30, 50, 70],
                'y': [-70, -50, -30, -10, 10, 30, 50, 70],
                'z': [-100, -75, -50, -25, 0, 25, 50, 75, 100]
                }
    for s in mykeys:
        print(f'Analyzing subject: {s}')
        for i in range(len(axes)):            
            axis = axes[i]
            plane = ''.join([k for k in poss_axes if k != axis])
            planeshape = myshape[plane]
            print(f'Axis: {axis}, corresponding plane: {plane}\nPlane size: {planeshape[0]} x {planeshape[1]} A')

            poss_axes = ['x', 'y', 'z']
            subj = f'{s}_{plane}'

            if densities:
                pd = outpath + '/densities/'
                try:
                    mkdir(pd)
                except OSError:
                    True 

                norm_mat = mats[subj]
                dname = pd + f'{subj}_density_{var}.png'
                densmap_plotter(norm_mat, dname, plane, cm, tick_dict)


            if len(contours) != 0:
                pc = outpath + '/contours/'
                try:
                    mkdir(pc)
                except OSError:
                    True 
                if s in targets:
                    #if s != 'prot_et_mem':
                    cname = pc + f'{subj}_TMhelices_contours_{var}.png'
                    contours_plotter(mats, contours, s, cname, plane, cm, tick_dict, colors = ['blue', 'green', 'purple', 'red', 'black'])



vars = ['wt', 'M654V'] #
axes = ['x', 'y', 'z']
contours = 'S4 S5 S6 S7 protein'.split()
targets = 'lipids potassium chloride prot_et_mem lipids'.split()
basepath = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/to_load/densmaps/' # str(sys.argv[1]) 
img_path = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/to_load/densmaps/images/' # str(sys.argv[2]) 

print('Starting')
for var in vars:
    path = f'{basepath}/{var}/'
    outpath = f'{img_path}/{var}/'
    try: 
        mkdir(outpath)
    except:
        True

    if var == 'wt':
        shape = {'xy':[140.0,  140.0], 'xz':[140.0,  200.0],'yz':[140.0,  200.0]}
    else:
        shape = {'xy':[140.0,  140.0], 'xz':[140.0,  202.0],'yz':[140.0,  202.0]}

    density_plotter(path, outpath, var, shape, targets, axes = axes, contours = contours, densities=True)
