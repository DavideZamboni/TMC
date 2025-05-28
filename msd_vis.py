import MDAnalysis as mda
from os import listdir, mkdir
import pandas as pd
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt

def smoothing(vals, times = None, smooth_window = 10):
    vals_std = []
    smt_vals = []
    start = 0
    if times is None:
        times = [i for i in range(len(vals))]
        
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
    
    return smt_vals, vals_std

def data_extraction(inpath, mylab = None):#, n_cols, labels):
    '''
    To get csvs from xvg files.
    Requires inpath
    '''
    print(inpath)
    labels = {}
    files = [f for f in listdir(inpath) if f.endswith('.xvg')]
    if mylab != None:
        k = [i for i in mylab.keys()]
        files = [f for f in listdir(inpath) if f.endswith('.xvg')]
        for ky in k:
            labels[ky] = {}
            actfiles = [f for f in files if ky in f]
            for f in actfiles:
                labels[ky][f] = []

    for key in k:
        files = [f for f in labels[key].keys()]
        for f in files:
            labels[ky][f] = ['x']
            #print(f, inpath + f)
            with open(inpath + f, 'r') as ff:
                a = ff.readlines()
                #print(a)
                if '&\n' in a:
                    print('multiple cols found')
                    for i in range(1, a.count('&\n') + 1):
                        labels[ky][f].append(f'y{i}')
                        print(labels[ky][f])
                else:
                    try:
                        inte = a[-1].split()
                        ncols = len(inte)
                        labels[ky][f] = ['x']
                        for i in range(1, ncols):
                            labels[ky][f].append(f'y{i}')
                    except IndexError as error:
                        break
    #print(labels)

    for key in k:
        files = [f for f in labels[key].keys()]
        for f in files:
            data = {}
            with open(inpath + f, 'r') as ff:
                inte = [l.split() for l in ff.readlines() if l[0] not in ['', '#', '@', '\n']]
                ncols = len(labels[ky][f])
                if ['&'] in inte:
                    i = 1
                    #print(ncols, labels[ky][f], inte.count(['&']))
                    for ll in inte:
                        if ll != ['&']:
                            if i == 1:
                                try:
                                    myx = float(ll[0])
                                    if 'x' not in data.keys():
                                        data['x'] = [myx]
                                    else:
                                        data['x'].append(myx)
                                except TypeError as error:
                                    #print('error: ', ll[0])
                                    if 'x' not in data.keys():
                                        data['x'] = [ll[0]]
                                    else:
                                        data['x'].append(ll[0])
                                try:
                                    myl = float(ll[1])
                                    if labels[ky][f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [myl]
                                    else:
                                        data[labels[ky][f][i]].append(myl)
                                except TypeError as error:
                                    #print('error: ', ll[1])
                                    if labels[ky][f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [ll[1]]
                                    else:
                                        data[labels[ky][f][i]].append(ll[1])
                            else:
                                try:
                                    myl = float(ll[1])
                                    if labels[ky][f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [myl]
                                    else:
                                        data[labels[f][i]].append(myl)
                                except TypeError as error:
                                    #print('error: ', ll[1])
                                    if labels[ky][f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [ll[1]]
                                    else:
                                        data[labels[ky][f][i]].append(ll[1])
                        elif ll == ['&']:
                            i += 1
                else:
                    for ll in inte:
                        if len(ll) != ncols:
                            print(f'warning, file {f} a different number of cols than expected, line: {"    ".join(ll)}')
                        if len(ll) == ncols:
                            i = 0
                            for l in ll:
                                try:
                                    myl = float(l)
                                    if labels[ky][f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [myl]
                                    else:
                                        data[labels[ky][f][i]].append(myl)
                                except TypeError as error:
                                    #print('error: ', l)
                                    if labels[f][i] not in data.keys():
                                        data[labels[ky][f][i]] = [myl]
                                    else:
                                        data[labels[f][i]].append(myl)
                                i += 1
            dat_df = pd.DataFrame.from_dict(data)
            #print(dat_df.head())
            if mylab != None:
                for ftype in mylab.keys():
                    if ftype in f:
                        if len(mylab[ftype]) == len(labels[ky][f]):
                            true_lab = {}
                            for i in range(len(labels[ky][f])):
                                true_lab[labels[ky][f][i]] = mylab[ftype][i]
                                #print(f'{f}: ', labels[ky][f][i], ' to ', mylab[ftype][i])
                            dat_df.rename(columns=true_lab, inplace=True)
                            #print(dat_df.head())
                            
                        else:
                            print('WARNING: ERROR IN LABELS')
            
            out = f'{inpath}resulting_dfs/'
            try:
                mkdir(out)
            except OSError as error:
                True
            
            dat_df.to_csv(f'{out}{f[:-3]}csv', sep = ',', index=False)
    return labels

def manipulation(inpath, labels, transform_var = None, smooth = None):
    '''
    transform_var: dict like {'rgyr':['Time (ns)', 'time (ps)', '/1000']} where at position 0: new variable, 1: old variable, 2: operation required to pass form old to new 
    smooth: dict like {'rgyr': []}
    '''
    out = f'{inpath}resulting_dfs/'
    for k in labels.keys():
        files = [f for f in labels[k].keys()]
        print(files)
        for f in files:
            if transform_var != None:
                keys = [i for i in transform_var.keys()]
                #print('The following measures will be transformed: ', keys)
                for k in keys:
                    if k in f:
                        if type(transform_var[k][0]) == list:
                            for new, old, operation in transform_var[k][:]:
                                df = pd.read_csv(f'{out}{f[:-3]}csv', index_col= None, header=0)
                                act_cols = [i for i in df.columns]
                                if old in act_cols:
                                    new_cols = [c if c != old else new for c in act_cols]
                                    if operation[0] == '+':
                                        df[new] = df[old] + float(operation[1:]) 
                                    elif operation[0] in ['*', 'x']:
                                        df[new] = df[old] * float(operation[1:]) 
                                    elif operation[0] == '-':
                                        df[new] = df[old] - float(operation[1:]) 
                                    elif operation[0] == '/':
                                        df[new] = df[old] / float(operation[1:]) 
                                    elif operation == 'count':
                                        df[new] = [i+1 for i in range(len(df[old]))]
                                        #print(df.head())
                                    df.to_csv(f'{out}{f[:-3]}csv', sep = ',', index=False, columns=new_cols)
                                else:
                                    print('No correspondence found, please check the spelling of the old variable: ', old)
                        else:
                            for new, old, operation in transform_var[k]:
                                df = pd.read_csv(f'{out}{f[:-3]}csv', index_col= None, header=0)
                                act_cols = [i for i in df.columns]
                                if old in act_cols:
                                    new_cols = [c if c != old else new for c in act_cols]
                                    if operation[0] == '+':
                                        df[new] = df[old] + float(operation[1:]) 

                                    elif operation[0] in ['*', 'x']:
                                        df[new] = df[old] * float(operation[1:]) 

                                    elif operation[0] == '-':
                                        df[new] = df[old] - float(operation[1:]) 

                                    elif operation[0] == '/':
                                        df[new] = df[old] / float(operation[1:]) 
                                    elif operation == 'count':
                                        df[new] = [i+1 for i in range(len(df[old]))]
                                        #print(df.head())
                                    df.to_csv(f'{out}{f[:-3]}csv', sep = ',', index=False, columns=new_cols)

                                else:
                                    print('No correspondence found, please check the spelling of the old variable: ', old)
            if smooth != None:
                keys = [i for i in smooth.keys()]
                for k in keys:
                    if k in f:
                        df = pd.read_csv(f'{out}{f[:-3]}csv', index_col= None, header=0)
                        act_cols = [i for i in df.columns]
                        try:
                            smooth_window = int(smooth[k][-1])
                        except TypeError as error:
                            print('Please insert an integer number as window size, now: ', smooth[k][-1])
                        if type(smooth[k][0]) == tuple or type(smooth[k][0]) == list:
                            for var in smooth[k][0]:
                                if var in act_cols and var + ' std' not in act_cols:
                                    vals = [i for i in df[var]]
                                    smt_vals, vals_std = smoothing(vals, smooth_window)
                                    df[var] = smt_vals
                                    df[var + ' std'] = vals_std
                                    #print(df.head())
                                    df.to_csv(f'{out}{f[:-3]}csv', sep = ',', index=False)
                                elif var in act_cols and var + ' std' in act_cols:
                                    print(f'In file {f}, the soothing operation has already been performed')
                                else:
                                    print('No correspondence found, please check the spelling of the variable: ', var)

                        elif type(smooth[k][0]) == str:
                            var = smooth[k][0]
                            #print(var)
                            if var in act_cols and var + ' std' not in act_cols:
                                vals = [i for i in df[var]]
                                smt_vals, vals_std = smoothing(vals, smooth_window)
                                df[var] = smt_vals
                                df[var + ' std'] = vals_std
                                #print(df.head())
                                df.to_csv(f'{out}{f[:-3]}csv', sep = ',', index=False)

                            elif var in act_cols and var + ' std' in act_cols:
                                print(f'In file {f}, the soothing operation has already been performed')
                            else:
                                print('No correspondence found, please check the spelling of the variable: ', var)

def concatenate(inpath, files, outname, update = None, add_col = None, skip=0):
    out = f'{inpath}resulting_dfs/'
    res_dic = {}
    i = 0
    for f in files:
        if f[:-3] + 'csv' not in listdir(out):
            print('Please perform data extraction before!')
        else:
            df = pd.read_csv(f'{out}{f[:-3]}csv', index_col= None, header=0)
            a = df.to_dict(orient='list')
            for k in a.keys():
                if k not in res_dic.keys():
                    res_dic[k] = a[k][skip:]
                    
                elif update != k:
                    for val in a[k][skip:]:
                        res_dic[k].append(val)
                elif update == k:
                    last = res_dic[k][-1]
                    step = mean([a[k][i + 1] - a[k][i] for i in range(len(a[k]) - 1)])
                    step = round(step, 3)
                    #print(last)
                    last -= step*skip
                    print(last)
                    #print(f'the identified step is {step} for the variable {k}')
                    print(res_dic[k][-10:])
                    for val in a[k][skip:]:
                        res_dic[k].append(round(val + last + step, 1))
            if add_col != None and add_col[0] not in res_dic.keys():
                res_dic[add_col[0]] = [add_col[1][i] for count in range(len(res_dic[k][skip:]))]
            elif add_col != None:
                for val in range(skip, len(a[k])):
                    res_dic[add_col[0]].append(add_col[1][i])
        i += 1
    enddf = pd.DataFrame.from_dict(res_dic)
    enddf.to_csv(f'{out}{outname}', sep = ',', index=False)

def get_pos_maps(df, struc, trj, ion):
    name = 'positions_w_lowD.csv'
    bp = '/'.join(df.split('/')[:-1])+'/'
    iname = ion[0]
    igro = ion[1]

    oname = 'lowD_gro_ids.txt'
    #if oname not in listdir(bp):
    u = mda.Universe(struc, trj)
    print(struc, ion)
    sdf = pd.read_csv(df, header=0, index_col=None)
    ion_ids = [i for i in u.select_atoms(f'name {igro}').ids]
    at_times = {'Time (ns)':[], 'GROMACS id':[]}
    for t, tdf in sdf.groupby(by='Time (ns)'):
        groids = [ion_ids[int(i)-1] for i in tdf['Atom id']]
        for groid in groids:
            at_times['Time (ns)'].append(int(t))
            at_times['GROMACS id'].append(int(groid))
    
    with open(bp + oname, 'w') as f:
        f.write('Time (ns),GROMACS id\n')
        for i in range(len(at_times['Time (ns)'])):
            f.write(f'{at_times["Time (ns)"][i]},{at_times["GROMACS id"][i]}\n')
    f.close()

    p = '/'.join(df.split('/')[:-1])+'/'
    #if name not in listdir(p):
    u = mda.Universe(struc, trj)
    sdf = pd.read_csv(df, header=0, index_col=None)
    ion_ids = [i for i in u.select_atoms(f'name {igro}').ids]
    start = 0
    lowD_pos = pd.DataFrame.from_dict({'Time (ns)':[], 'x':[], 'y':[], 'z':[], 'loc':[]})
    for t, tdf in sdf.groupby(by='Time (ns)'):
        myids = f'name {igro} and (id ' + ' or id '.join([f'{ion_ids[int(i)-1]}' for i in tdf['Atom id']]) + ')'
        print(myids)
        myframes = [ts.frame for ts in u.trajectory if start*1000 < ts.time <= t*1000]
        print(myframes)
        ats = u.select_atoms(myids)
        for ts in u.trajectory[myframes]:
            ps = ats.positions
            ids = ats.ids
            tdf = pd.DataFrame.from_dict({'Time (ns)':[ts.time/1000 for i in range(len(ps[:, 0]))], 'x':ps[:, 0], 'y':ps[:, 1], 'z':ps[:, 2], 'loc':[iname+str(ids[i]) for i in range(len(ps[:, 0]))]})
            lowD_pos = pd.concat([lowD_pos, tdf])
        start = t
    lowD_pos.to_csv(p + name, header = True, index=False)

def get_data(path, names):
    '''
    
    '''
    ion = names[2]
    observables = {'diff_mol':[f'Mean Diffusion coefficient for {ion}', {'D':['Diffusion Coefficient ' + r'$\left(10^{-5}  \frac{cm²}{s}\right)$']}, 'Atom id', 'D ' + r'$\left(10^{-5}  \frac{cm²}{s}\right)$'],
                   'msdout':[f'Mean Squared Displacement (nm²) for {ion}', {'trace':'Tot', 'xx':'x-axis', 'yy':'y-axis', 'zz':'z-axis', 'yx':'xy-plane', 'zx':'xz-plane', 'zy':'yz-plane'}, 'Time (ns)', 'MSD(t)']}      

    labels = {'diff_mol':['Atom id', 'Diffusion Coefficient (cm²/s)'],
              'msdout':['time (ps)', 'trace', 'xx', 'yy', 'zz', 'yx', 'zx', 'zy']}
    mylabs = data_extraction(path, labels)
    trans_var =  {'msdout':[['Time (ns)', 'time (ps)', '/1000']], 'diff_mol_':[['Atom id', 'Atom id', 'count']]}
    manipulation(path, mylabs, transform_var=trans_var)
    for o in observables.keys():
        try:
            pp = path + '/resulting_dfs/'
            mkdir(pp)
        except OSError:
            True
        concat_files = [f for f in listdir(path + '/resulting_dfs/') if o in f and not 'merged' in f]
        outname = f'{o}_merged.csv'
        if f'{o}_merged.csv' not in listdir(path + '/resulting_dfs/'):
            if o == 'diff_mol':
                times = [int(f.split('_')[-1][4:-4])/1000 for f in concat_files if o in f and '_merged.csv' not in f]
                times.sort()
                print(times)
                ac = ['Time (ns)', times]
                concatenate(path, concat_files, outname, add_col = ac)
                continue
            elif o == 'msdout':
                concatenate(path, concat_files, outname, update= 'Time (ns)', skip=1)
                continue

def msd_plots(path, pdir, threshold = 0.05, diff_threshold = None):
    inname = f'{path}/resulting_dfs/diff_mol_merged.csv'
    diff_df = pd.read_csv(f'{inname}', header = 0, index_col=False)

    print(inname)
    p = '/'.join(inname.split('/')[:-1]) + '/'

    count, division = np.histogram(diff_df['Diffusion Coefficient (cm²/s)'], bins = 1000, density=False)
    sums = np.cumsum(count)
    densities = [i/sums[-1] for i in count]
    cumulative_sum = np.cumsum(densities)
    last_bin_index = np.searchsorted(cumulative_sum, threshold, side='right') - 1
    if diff_threshold is None:
        diff_threshold = division[last_bin_index]

    small=8
    medium=16
    large=18
    xl=22
    plt.rc('font',family="serif")
    plt.rc('axes', labelsize=small)
    plt.rc('xtick', labelsize=small)
    plt.rc('ytick', labelsize=small)
    plt.rc('legend', fontsize=small)
    plt.hist(diff_df['Diffusion Coefficient (cm²/s)'], bins = 1000, density=True)
    # Plot a vertical line at the threshold
    plt.axvline(x=diff_threshold, color='r', linestyle='--', label=f'Threshold ({diff_threshold})')
    print(f"Diffusion coefficient threshold to identify the {threshold*100}% slow-diffusing ions: {diff_threshold}")
    plt.ylabel('Frequency')
    plt.xlabel('D '+ r'$\left(\frac{cm²}{s}\right)$')
    plt.xlim(0, 30)
    plt.ylim(0, 0.30)
    plt.legend()
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f'{pdir}Diffusion_coeff_distribution_{threshold}.png', dpi = 600)
    plt.close()
    bins = [division[i] for i in range(len(division)-1)]
    print(len(bins))
    myd = {'Bins':bins, 'Frequency':count, 'Density':densities, 'Cumulative density':cumulative_sum}
    
    for k in myd.keys():
        print(k, len(myd[k]))
    data = pd.DataFrame.from_dict(myd)
    data.to_csv(p + 'diffusion_distribution.csv', header=True, index=False)
    
    sdf = pd.DataFrame.from_dict({'Time (ns)':[], 'Atom id':[], 'Diffusion Coefficient (cm²/s)':[]})
    dtypes = {'Time (ns)': float, 'Atom id': int, 'Diffusion Coefficient (cm²/s)': float}
    # Assign data types to columns
    sdf = sdf.astype(dtypes)

    for t, tdf in diff_df.groupby(by = 'Time (ns)'):
        minif = tdf[tdf['Diffusion Coefficient (cm²/s)'] <= diff_threshold]
        sdf = pd.concat([sdf, minif])
    f = inname.split('/')[-1].replace('merged', 'slowdiff_ions')
    sdf.to_csv(p + f, header=True, index=False)

    return diff_threshold


bp = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/to_load/msd/'
pdd = f'{bp}images/'
vars = ['wt', 'M654V']  # 
locs = ['potassium', 'chloride'] 
thresholds = {}
#### TO ANALYZE IONS DISTRIBUTION OF THE SLOW-DIFFUSING IONS
for var in vars:
    thresholds[var] = {}
    try:
        mkdir(f'{pdd}{var}')
    except OSError:
        True
    for loc in locs:
        pdir = f'{pdd}{var}/{loc}/'
        try:
            mkdir(pdir)
        except OSError:
            True

        if loc == 'chloride':
            names = [var, loc, r'$Cl^{-}$']
            ion = ['chloride', 'CLA', r'$Cl^{-}$']
        elif loc == 'potassium':
            names = [var, loc, r'$K^{+}$']
            ion = ['potassium', 'POT', r'$K^{+}$']
        p = f'{bp}/{var}/{loc}/'
        trjpath = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/trj_data/{var}/xtcs/'
        get_data(p, names)
        thresholds[var][loc] = msd_plots(p, pdir, threshold = 0.05)

        df = f'{p}resulting_dfs/diff_mol_slowdiff_ions.csv'
        trj = f'{trjpath}TMC1_{var}_noPBC.xtc'
        stru = f'{trjpath}TMC1_{var}_noPBC.gro'
        u = mda.Universe(stru, trj)

        get_pos_maps(df, stru, trj, ion)

