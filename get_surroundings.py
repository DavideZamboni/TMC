import MDAnalysis as mda
import numpy as np
import pandas as pd
from statistics import mean, stdev
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_surroundings(struc, trj, resid, out, radius = 5, monomers = ['A', 'B'], monomer_residues = 760, exclude = 2):
    '''
    
    Selects all atoms within a spherical zone centered in the center of geometry (COG) of the side chain of a given residue resid.
    Writes a .csv with columns: 'Monomer', 'Time (ns)', 'Resname', 'Resid'
    '''
    u = mda.Universe(struc, trj)
    results = {'Monomer':[], 'Time (ns)':[], 'Resname':[], 'Resid':[]}
    for i in range(len(monomers)):
        idx = resid + i * monomer_residues

        surrounding = u.select_atoms(f'sphzone {radius} (resid {idx} and not backbone)', updating=True)
        for ts in u.trajectory:
            just_added = []
            for res in surrounding.residues:
                if abs(res.resid - idx) > exclude and res.resid not in just_added:
                    results['Monomer'].append(monomers[i])
                    results['Time (ns)'].append(ts.time/1000)
                    results['Resname'].append(res.resname)
                    results['Resid'].append(res.resid)
                    just_added.append(res.resid)
    df = pd.DataFrame.from_dict(results)
    #print(df.head(50))
    df.to_csv(out, header = True, index=False)

def compare_w_minrad(out_surroundings, hole_features, out, monomers = ['A', 'B']):
    '''
    Compares the data from get_surroundings with minrad data from HOLE analysis.
    In particular, counts the times that any residue found in the surroundings of the amalyzed zone
    is found at the point of minimum radius of the pore in such timestep.

    Writes a .txt file with:
    - header: a section exportable as a csv with columns Resname,Resid,Time (ns)
    - At the bottom, the proportion of resiudes close to mutation and chocke point and the total number of timesteps checked.
    '''
    surroundings = pd.read_csv(out_surroundings, header=0, index_col=False)
    for i  in range(len(monomers)):
        surs = surroundings[surroundings['Monomer'] == monomers[i]]
        features = pd.read_csv(hole_features[i], header=0, index_col=False)
        checked = 0
        results = []
        
        results.append(f"Resname,Resid,Time (ns)")

        for t, tfeats in features.groupby(by = 'Time (ns)'):
            #print(t)
            minrad = np.min(tfeats['Radius $\AA$'])
            #print(t, minrad)
            ms = surs[surs['Time (ns)'] == t]
            minress = tfeats[tfeats['Radius $\AA$'] == minrad]['Resid']
            #print(minress)
            myress = ms[ms['Resid'].isin(minress)][['Resid', 'Resname']]
            #print(myress.head())
            checked += 1
            if myress.shape[0] > 0:
                for ix, r in myress.iterrows():
                    results.append(f"{r['Resname']},{r['Resid']},{t}")
        results.append(f'Proportion of Resiudes close to mutation and chocke point: {len(results)/checked}')
        results.append(f'Tot checked: {checked}')
        with open(out[i], 'w') as f:
            f.write('\n'.join(results))

def get_helix_z_distances(struc, trj, resid, out, myhelix = 'S9', z_extension = 2.5, monomers = ['A', 'B'], monomer_residues = 760, extend = 5, 
                        tm_helices = {'S1':(180,220), 'S2':(280,320), 'S3':(350,380),
                        'S4':(402,427), 'S5':(435,465), 'S6':(517,552),
                        'S7':(574,592), 'S8':(595,620), 'S9':(635,655),
                        'S10':(700,730)}
                        ):
    
    '''
    Selects all atoms within a spherical zone centered in the center of geometry (COG) of a given residue resid.
    Writes a .csv with columns: 'Monomer', 'Time (ns)', 'Resname', 'Resid'
    '''
    
    distances = {'Monomer':[], 'Time (ns)':[]}
    for k in tm_helices.keys():
        if k != myhelix:
            distances[k] = []

    u = mda.Universe(struc, trj)
    for i in range(len(monomers)):
        idx = resid + i * monomer_residues
        myres = u.select_atoms(f'protein and resid {idx}', updating=True)
        for ts in u.trajectory:
            #res_cog = myres.center_of_geometry()
            res_com = myres.center_of_mass()
            distances['Time (ns)'].append(ts.time/1000)
            distances['Monomer'].append(monomers[i])
            for h in tm_helices.keys():
                if h != myhelix:
                    helix_sel = u.select_atoms(f'''protein and (resid {tm_helices[h][0] + i * monomer_residues + extend}:{tm_helices[h][1] + i * monomer_residues + extend}) 
                                               and (prop z >= {res_com[2] - z_extension}) and (prop z <= {res_com[2] + z_extension})''')
                    m = 2
                    while len(helix_sel.atoms) < 1:
                        #print(f'extending z_extension to {m * z_extension} for helix {h} in monomer {monomers[i]} at time {ts.time/1000} ns')
                        helix_sel = u.select_atoms(f'''protein and (resid {tm_helices[h][0] + i * monomer_residues + extend}:{tm_helices[h][1] + i * monomer_residues + extend})
                                                   and (prop z >= {res_com[2] - m * z_extension}) and (prop z <= {res_com[2] + m * z_extension})''')
                        m += 1
                    hcog = helix_sel.center_of_geometry()
                    hcom = helix_sel.center_of_mass()
                    com_dist = np.linalg.norm(res_com - hcom)
                    #cog_dist = np.linalg.norm(res_cog - hcog)
                    distances[h].append(com_dist)

    df = pd.DataFrame.from_dict(distances)
    df.to_csv(out, header = True, index=False)

def get_helix_cog_distances(struc, trj, resid, out, myhelix = 'S9', extension = 2.5, monomers = ['A', 'B'], monomer_residues = 760, extend = 5, 
                        tm_helices = {'S1':(180,220), 'S2':(280,320), 'S3':(350,380),
                        'S4':(402,427), 'S5':(435,465), 'S6':(517,552),
                        'S7':(574,592), 'S8':(595,620), 'S9':(635,655),
                        'S10':(700,730)}, backbone = True
                        ):
    
    '''
    Selects all atoms within a spherical zone centered in the center of geometry (COG) of a given residue resid.
    Writes a .csv with columns: 'Monomer', 'Time (ns)', 'Resname', 'Resid'
    '''
    
    distances = {'Monomer':[], 'Time (ns)':[]}
    for k in tm_helices.keys():
        if k != myhelix:
            distances[k] = []

    u = mda.Universe(struc, trj)
    for i in range(len(monomers)):
        idx = resid + i * monomer_residues
        if backbone:
            myres = u.select_atoms(f'protein and resid {idx}', updating=True)
        else:
            myres = u.select_atoms(f'protein and resid {idx} and not backbone', updating=True)

        for ts in u.trajectory:
            res_cog = myres.center_of_geometry()
            res_com = myres.center_of_mass()
            distances['Time (ns)'].append(ts.time/1000)
            distances['Monomer'].append(monomers[i])
            for h in tm_helices.keys():
                if h != myhelix:
                    helix_sel = u.select_atoms(f'''protein and (resid {tm_helices[h][0] + i * monomer_residues + extend}:{tm_helices[h][1] + i * monomer_residues + extend}) 
                                               and (sphzone {extension} (protein and resid {idx}))''')
                    m = 2
                    while len(helix_sel.atoms) < 1:
                        #print(f'extending tadis to {m * extension} for helix {h} in monomer {monomers[i]} at time {ts.time/1000} ns')
                        helix_sel = u.select_atoms(f'''protein and (resid {tm_helices[h][0] + i * monomer_residues + extend}:{tm_helices[h][1] + i * monomer_residues + extend})
                                                       and (sphzone {m * extension} (protein and resid {idx}))''')
                        m += 1
                    rescogs = [r.atoms.center_of_geometry() for r in helix_sel.residues]
                    cog_dist = np.min([np.linalg.norm(res_cog - hcog) for hcog in rescogs])
                    distances[h].append(cog_dist)

    df = pd.DataFrame.from_dict(distances)
    df.to_csv(out, header = True, index=False)

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

def compare_dists_byhelix(outs, outplot_path, vars = ['wt', 'M654V'], tm_helices=['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10'], ref = 'S9', smooth_window = 10):
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

    for h in tm_helices:
        if h != ref:
            color_count = 0
            for i in range(len(outs)):
                heldf = pd.read_csv(outs[i], header=0, index_col=False, usecols=['Time (ns)', 'Monomer', h])

                for mon, df in heldf.groupby(by = 'Monomer'):
                    df = df.sort_values('Time (ns)')

                    df['Smooth Distance '+ h + r'$\AA$'], df['Smooth Distance '+ h + ' std'] = time_smoothing(df[h].to_list(), df['Time (ns)'].to_list(), smooth_window)
                    plt.plot(df['Time (ns)'], df['Smooth Distance '+ h + r'$\AA$'], label = f'{mon} {vars[i]}', color=colors[color_count])
                    plt.fill_between(df['Time (ns)'], df['Smooth Distance '+ h + r'$\AA$']-df['Smooth Distance '+ h + ' std'], df['Smooth Distance '+ h + r'$\AA$']+df['Smooth Distance '+ h + ' std'], color=colors[color_count], alpha=0.5)
                    color_count += 1
            plt.xlabel('Time (ns)')
            plt.ylabel('Distance ('+ r'$\AA$' + ')')
            leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
            leg.get_frame().set_linewidth(0.0)  
            plt.xlim(0, 1000)
            #plt.ylim(0, 2)
            plt.xticks(range(0, 1001, 200))
            #plt.yticks(np.linspace(0, 2, 9))
            plt.tight_layout()
            plt.savefig(outplot_path + f'{h}_cog_distance_vs_time.png', dpi=600)
            plt.close()

def interactors_histograms(outs, vars, outplot_path, monomer_residues = 760, min_freq = 5):
    colors = ["dodgerblue","red", "forestgreen", "gold", "darkviolet", "black"]
    small=12
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

    results = {}
    monomers = set()
    # Compute and plot % frequency of appearance of resname resindex (on the number of analyzed frames) per variant per monomer
    # Compute the results for each variant and plot % frequency of appearance of resname resindex (on the number of analyzed frames)
    for i in range(len(vars)):
        results[vars[i]] = {}
        df = pd.read_csv(outs[i], header=0, index_col=False)
        nmons = len(np.unique(df['Monomer']))
        #pervar_ntimes = len(np.unique(df['Time (ns)']))
        for m, mdf in df.groupby(by = 'Monomer'):
            monomers.add(m)
            ntimes = len(np.unique(mdf['Time (ns)']))
            resfreq = []
            for ri, ridf in mdf.groupby(by = 'Resid'):
                rname = np.unique(ridf['Resname'])
                #print(rname)
                if len(rname) == 1:
                    if ri <= monomer_residues or rname[0] == 'POPC':
                        myid = ri
                    else:
                        myid = ri-monomer_residues
                    res = f'{rname[0]} {myid}'
                    resfreq.append([res, 100 * ridf.shape[0]/ntimes])
                else:
                    print(f'Warning: more resnames associated to the residue index {ri}: {rname}')
            resfreq.sort(key=lambda x: x[1], reverse = True)
            myfreq = [i for i in resfreq if i[1] >= min_freq ]
            results[vars[i]][m] = resfreq
            plt.title(f'Monomer {m}, TMC1 {vars[i]}')
            plt.bar([i[0] for i in myfreq], [i[1] for i in myfreq])
            plt.xticks(rotation = 90)
            plt.ylabel('% frequency')
            plt.tight_layout()
            plt.savefig(outplot_path + f'interactors_within_5_{vars[i]}_{m}.png', dpi=600)
            plt.close()

    per_variant_results = {}

    for v in results.keys():
        already_found = []
        per_variant_results[v] = []
        nmons = len(monomers)
        for m in results[v].keys():
            for [res, freq] in results[v][m]:
                #print(res, freq)
                if res in already_found:
                    per_variant_results[v][already_found.index(res)][1] += freq/nmons
                    #print(per_variant_results[v][already_found.index(res)][1])

                else: 
                    already_found.append(res)
                    per_variant_results[v].append([res, freq/nmons])


        per_variant_results[v].sort(key=lambda x: x[1], reverse = True)
        freqs = [i for i in per_variant_results[v] if i[1] >= min_freq]
        plt.title(f'TMC1 {v}')
        plt.bar([i[0] for i in freqs], [i[1] for i in freqs])
        plt.xticks(rotation = 90)
        plt.ylabel('% frequency')
        plt.tight_layout()
        plt.savefig(outplot_path + f'interactors_within_5_{v}_tot.png', dpi=600)
        plt.close()

    res = {'Variant':[], 'Residue':[], 'Freq. (%)':[]}
    
    for v in results.keys():
        #print(v, per_variant_results[v])
        for (resi, freq) in per_variant_results[v]:
            #print((res, freq))
            res['Variant'].append(v)
            res['Residue'].append(resi)
            res['Freq. (%)'].append(freq)
    
    df = pd.DataFrame.from_dict(res)
    df.to_csv(outplot_path + 'Pervariant_interactors_freqs.csv', header = True, index=False)

    res = {'Variant':[], 'Monomer':[], 'Residue':[], 'Freq. (%)':[]}
    
    for v in results.keys():
        for m in results[v].keys():
            for (resi, freq) in results[v][m]:
                #print((res, freq))
                res['Variant'].append(v)
                res['Monomer'].append(m)
                res['Residue'].append(resi)
                res['Freq. (%)'].append(freq)
    
    df = pd.DataFrame.from_dict(res)
    df.to_csv(outplot_path + 'Per_monomer_interactors_freqs.csv', header = True, index=False)


    # Plot differences in % frequency, vars[0] - vars[1], vars[0] - vars[2], vars[0] - vars[3], vars[1] - vars[2], ... per var and per var per monomer
    # Find unique residues (residues that appears only in one variant, regardless of the monomer)

    ## Per monomer differences
    unique = {}
    for v in vars:
        unique[v] = []

    for m in monomers:
        differences = []
        
        for i in range(len(vars)):
            for j in range(i+1, len(vars)):
                results_i = results[vars[i]][m]
                results_j = results[vars[j]][m]
                residues_i = [k[0] for k in results_i]
                residues_j = [k[0] for k in results_j]
                for res in residues_i:
                    if res in residues_j:
                        freq_i = results_i[residues_i.index(res)][1]
                        freq_j = results_j[residues_j.index(res)][1]
                        if abs(freq_i - freq_j) >= min_freq:
                            differences.append((res, freq_i - freq_j))
                    else:
                        unique[vars[i]].append(res)
                        freq_i = results_i[residues_i.index(res)][1]
                        differences.append((res, freq_i))
                for res in residues_j:
                    if res not in residues_i:
                        unique[vars[j]].append(res)

                        freq_j = results_j[residues_j.index(res)][1]
                        differences.append((res, -freq_j))

                differences.sort(key=lambda x: x[1], reverse = True)
                #print(m, differences)
                #plt.title(f'Monomer {m} differences')
                colors = ['dodgerblue' if freq[1] > 0 else 'red' for freq in differences]

                plt.bar([u[0] for u in differences], [u[1] for u in differences], color=colors)
                plt.xticks(rotation = 90)
                plt.ylabel('% frequency')
                plt.axhline(y = 0, color='black', linestyle='--', zorder = 1)
                xlims = plt.xlim()
                labels = [[round(xlims[1]/2) + 1, 80, f'Present in {vars[i]}', 'dodgerblue'], 
                          [1, -80, f'Present in {vars[j]}', 'red']]

                for label in labels:
                    plt.text(label[0], label[1], label[2], fontsize=small, color=label[3])

                plt.ylim((-100, 100))
                plt.tight_layout()
                #plt.show()
                plt.savefig(outplot_path + f'differences_monomer_{m}.png', dpi=600)
                plt.close()

    ## Per variant differences
    differences = []
    for i in range(len(vars)):
        for j in range(i+1, len(vars)):
            results_i = per_variant_results[vars[i]]
            results_j = per_variant_results[vars[j]]
            residues_i = [k[0] for k in results_i]
            residues_j = [k[0] for k in results_j]
            for res in residues_i:
                if res in residues_j:
                    freq_i = results_i[residues_i.index(res)][1]
                    freq_j = results_j[residues_j.index(res)][1]
                    if abs(freq_i - freq_j) >= min_freq:
                        differences.append((res, freq_i - freq_j))
                else:
                    freq_i = results_i[residues_i.index(res)][1]
                    differences.append((res, freq_i))
            for res in residues_j:
                if res not in residues_i:
                    freq_j = results_j[residues_j.index(res)][1]
                    differences.append((res, -freq_j))
            
            #print(differences)
            differences.sort(key=lambda x: x[1], reverse = True)
            mydiffs = [k for k in differences if abs(k[1]) >= min_freq]
            #print(vars[i], vars[j], differences)
            colors = ['dodgerblue' if freq[1] > 0 else 'red' for freq in mydiffs]
            plt.bar([u[0] for u in mydiffs], [u[1] for u in mydiffs], color=colors)
            plt.xticks(rotation = 90)
            plt.ylabel('% frequency')
            plt.axhline(y = 0, color='black', linestyle='--', zorder = 1)
            xlims = plt.xlim()
            labels = [[round(xlims[1]/2) + 1, 80, f'Present in {vars[i]}', 'dodgerblue'], 
                      [1, -80, f'Present in {vars[j]}', 'red']]
            for label in labels:
                plt.text(label[0], label[1], label[2], fontsize=small, color=label[3])

            plt.ylim((-100, 100))
            plt.tight_layout()
            #plt.show()
            plt.savefig(outplot_path + f'differences_per_variant_{vars[i]}minus{vars[j]}.png', dpi=600)
            plt.close()

            d = {'Variant':[], 'Residue':[], 'Freq (%)':[]}
            for (r, f) in differences:
                if f > 0:
                    d['Variant'].append(vars[i])
                else:
                    d['Variant'].append(vars[j])
                d['Residue'].append(r)
                d['Freq (%)'].append(f)
                df = pd.DataFrame.from_dict(d)
                df.to_csv(outplot_path + f'{vars[i]}-{vars[j]}_interactors_diffs.csv', header = True, index=False)

    for i in range(len(vars)):
        for j in range(i+1, len(vars)):
            res_i = unique[vars[i]]
            res_j = unique[vars[j]]
            for r in res_i:
                if r in res_j:
                    unique[vars[i]].remove(r)
                    unique[vars[j]].remove(r)
            for r in res_j:
                if r in res_i:
                    unique[vars[i]].remove(r)
                    unique[vars[j]].remove(r)
            unique[vars[i]] = set(unique[vars[i]])
            unique[vars[j]] = set(unique[vars[j]])
    print('unique', unique)

def get_aver_position(struc, trj, resid, out, monomers = ['A', 'B'], monomer_residues = 760, backbone = True):
    
    '''
    Selects all atoms within a spherical zone centered in the center of geometry (COG) of a given residue resid.
    Writes a .csv with columns: 'Monomer', 'Time (ns)', 'Resname', 'Resid'
    '''
    
    positions = {'Monomer':[], 'Time (ns)':[], 'x':[], 'y':[], 'z':[]}

    u = mda.Universe(struc, trj)
    for i in range(len(monomers)):
        idx = resid + i * monomer_residues
        if backbone:
            myres = u.select_atoms(f'protein and resid {idx}', updating=True)
        else:
            myres = u.select_atoms(f'protein and resid {idx} and not backbone', updating=True)

        for ts in u.trajectory:
            res_cog = myres.center_of_geometry()
            positions['Time (ns)'].append(ts.time/1000)
            positions['Monomer'].append(monomers[i])
            positions['x'].append(res_cog[0])
            positions['y'].append(res_cog[1])
            positions['z'].append(res_cog[2])

    res = pd.DataFrame.from_dict(positions)
    avgs = {'Monomer':[], 'Time (ns)':[], 'x':[], 'y':[], 'z':[]}
    for m, mdf in res.groupby(by = 'Monomer'):
        avx = np.mean(mdf['x'].values)
        avy = np.mean(mdf['y'].values)
        avz = np.mean(mdf['z'].values)
        avgs['Time (ns)'].append('trj')
        avgs['Monomer'].append(f'Average {m}')
        avgs['x'].append(avx)
        avgs['y'].append(avy)
        avgs['z'].append(avz)

        stdx = np.std(mdf['x'].values)
        stdy = np.std(mdf['y'].values)
        stdz = np.std(mdf['z'].values)
        avgs['Time (ns)'].append('trj')
        avgs['Monomer'].append(f'Stdev {m}')
        avgs['x'].append(stdx)
        avgs['y'].append(stdy)
        avgs['z'].append(stdz)
    
    res = pd.concat([res, pd.DataFrame.from_dict(avgs)], ignore_index=True)
    res.to_csv(out, header = True, index=False)

def assign_regions(data, tm_helices = {'S1':(180,220), 'S2':(280,320), 'S3':(350,380),
                        'S4':(402,427), 'S5':(435,465), 'S6':(517,552),
                        'S7':(574,592), 'S8':(595,620), 'S9':(635,655),
                        'S10':(700,730)}, maxid = 760
                        ):
    hels = [k for k in tm_helices.keys()]
    hels.sort(key = lambda x: tm_helices[x][0])
    #print(hels)

    d = pd.read_csv(data, header=0, index_col=False)
    for r, rdf in d.groupby(by='Residue'):
        idx = int(r.split(' ')[1])
        for i in range(len(hels)):
            if tm_helices[hels[i]][0] <= idx <= tm_helices[hels[i]][1]:
                d.loc[rdf.index, 'Region'] = hels[i]
            elif i == 0 and idx < tm_helices[hels[i]][0]:
                d.loc[rdf.index, 'Region'] = 'beforeS1'
            elif i == len(hels) - 1 and idx > tm_helices[hels[i]][1]:
                d.loc[rdf.index, 'Region'] = 'afterS10'
            elif tm_helices[hels[i]][1] < idx < tm_helices[hels[i + 1]][0]:
                d.loc[rdf.index, 'Region'] = f'S{hels[i][1:]}/{hels[i+1][1:]}'
    d.to_csv(data, header = True, index=False)


def xy_radgyr(atomgroup, masses, total_mass=None):
    # Coordinates of the atoms (changes per frame)
    coordinates = atomgroup.positions  # shape: (n_atoms, 3)
    
    # Center of mass (already mass-weighted in MDAnalysis)
    center_of_mass = atomgroup.center_of_mass()
    
    # Squared distance from COM (for x, y only)
    ri_sq_xy = (coordinates[:, :2] - center_of_mass[:2])**2  # shape: (n_atoms, 2)
    sq_xy = np.sum(ri_sq_xy, axis=1)  # sum x² + y² for each atom
        
    # Mass-weighted Rg² in the xy-plane
    rog_sq = np.sum(masses * sq_xy) / total_mass
    
    # Return Rg (sqrt of Rg²)
    return np.sqrt(rog_sq)

def get_xy_gyration(struc, trj, out, monomers=['A', 'B'], maxresid = 760, tm_helices = {'S1':(182,212), 'S2':(274,303), 'S3':(347,392),
                          'S4':(402,427), 'S5':(435,462), 'S6':(517,552,
                          'S7':(574,589), 'S8':(595,620), 'S9':(635,665),
                          'S10':(700,730)}):

    u = mda.Universe(struc, trj)
    d = {'Time (ns)':[], 'Monomer':[], 'xy-rgyr, (' + r'$\AA$' + ')':[]}
    # per monomer:
    whole = 'protein and ('
    for i in range(len(monomers)):
        sel = 'protein and ('
        for k in tm_helices.keys():
            sel += f'resid {tm_helices[k][0] + i * maxresid}:{tm_helices[k][1] + i * maxresid} or '
            whole += f'resid {tm_helices[k][0] + i * maxresid}:{tm_helices[k][1] + i * maxresid} or '
             
        sel = sel[:-4] + ')'
        ats = u.select_atoms(sel)
        masses = ats.masses
        total_mass = np.sum(ats.masses)
        m = monomers[i]
        for ts in u.trajectory:
            xyrgyr = xy_radgyr(ats, masses, total_mass=total_mass)
            d['Time (ns)'].append(ts.time / 1000)
            d['Monomer'].append(m)
            d['xy-rgyr, (' + r'$\AA$' + ')'].append(xyrgyr)
    
    
    whole = whole[:-4] + ')'
    ats = u.select_atoms(whole)
    masses = ats.masses
    total_mass = np.sum(ats.masses)
    m = monomers[i]
    for ts in u.trajectory:
        xyrgyr = xy_radgyr(ats, masses, total_mass=total_mass)
        d['Time (ns)'].append(ts.time / 1000)
        d['Monomer'].append('Total')
        d['xy-rgyr, (' + r'$\AA$' + ')'].append(xyrgyr)

    df = pd.DataFrame.from_dict(d)
    df.to_csv(out, header = True, index=False)

def compare_rgyrs(outs, vars, outplot_path, smooth_window = 10, monomers = ['A', 'B']):
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

    color_count = 0
    nonmons = set()
    for i in range(len(outs)):
        heldf = pd.read_csv(outs[i], header=0, index_col=False)
        for mon, df in heldf.groupby(by = 'Monomer'):
            if mon in monomers:
                df = df.sort_values('Time (ns)')
                df['xy-rgyr, (' + r'$\AA$' + ')'], df['Smooth Distance xy-rgyr std'] = time_smoothing(df['xy-rgyr, (' + r'$\AA$' + ')'].to_list(), df['Time (ns)'].to_list(), smooth_window)
                plt.plot(df['Time (ns)'], df['xy-rgyr, (' + r'$\AA$' + ')'], label = f'{mon} {vars[i]}', color=colors[color_count])
                plt.fill_between(df['Time (ns)'], df['xy-rgyr, (' + r'$\AA$' + ')']-df['Smooth Distance xy-rgyr std'], df['xy-rgyr, (' + r'$\AA$' + ')']+df['Smooth Distance xy-rgyr std'], color=colors[color_count], alpha=0.5)
                color_count += 1
            else:
                nonmons.add(mon)
    plt.xlabel('Time (ns)')
    plt.ylabel('xy-'+ r'$R_{gyr}$' + ', (' + r'$\AA$' + ')')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(0, 1000)
    #plt.ylim(0, 2)
    plt.xticks(range(0, 1001, 200))
    #plt.yticks(np.linspace(0, 2, 9))
    plt.tight_layout()
    plt.savefig(outplot_path + 'permon_comparison_xyrgyr_vs_time.png', dpi=600)
    plt.close()

    color_count = 0
    for i in range(len(outs)):
        heldf = pd.read_csv(outs[i], header=0, index_col=False)
        for mon, df in heldf.groupby(by = 'Monomer'):
            if mon in nonmons:
                df = df.sort_values('Time (ns)')
                df['xy-rgyr, (' + r'$\AA$' + ')'], df['Smooth Distance xy-rgyr std'] = time_smoothing(df['xy-rgyr, (' + r'$\AA$' + ')'].to_list(), df['Time (ns)'].to_list(), smooth_window)
                plt.plot(df['Time (ns)'], df['xy-rgyr, (' + r'$\AA$' + ')'], label = f'{mon} {vars[i]}', color=colors[color_count])
                plt.fill_between(df['Time (ns)'], df['xy-rgyr, (' + r'$\AA$' + ')']-df['Smooth Distance xy-rgyr std'], df['xy-rgyr, (' + r'$\AA$' + ')']+df['Smooth Distance xy-rgyr std'], color=colors[color_count], alpha=0.5)
                color_count += 1
    plt.xlabel('Time (ns)')
    plt.ylabel('xy-' + r'$R_{gyr}$' + ', (' + r'$\AA$' + ')')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(0, 1000)
    #plt.ylim(0, 2)
    plt.xticks(range(0, 1001, 200))
    #plt.yticks(np.linspace(0, 2, 9))
    plt.tight_layout()
    plt.savefig(outplot_path + 'Total_comparison_xyrgyr_vs_time.png', dpi=600)
    plt.close()

def get_avgdists(outs, vars, outplot_path, smooth_window = 10):
    ## Average distance from S9 to all other helices at the mutation point
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

    helices = [f'S{i}' for i in range(1, 11) if i != 9]
    tot_avgs = {'Variant':[], 'Time (ns)':[], 'Avg dist':[]}
    for i in range(len(outs)):
        heldf = pd.read_csv(outs[i], header=0, index_col=False)
        heldf = heldf.sort_values('Time (ns)')

        for t, tdf in heldf.groupby(by = 'Time (ns)'):
            tot_avgs['Time (ns)'].append(t)
            tot_avgs['Variant'].append(vars[i])
            alld = 0
            for ix, row in tdf.iterrows():
                mymean = row[helices].mean(axis = 0)
                heldf.loc[ix, 'Avg dist'] = mymean
                alld += mymean
            tot_avgs['Avg dist'].append(alld/tdf.shape[0])

        heldf.to_csv(outs[i], header = True, index=False)

    '''color_count = 0
    for i in range(len(outs)):
        heldf = pd.read_csv(outs[i], header=0, index_col=False, usecols=['Monomer', 'Time (ns)', 'Avg dist'])
        heldf = heldf.sort_values('Time (ns)')
        for mon, df in heldf.groupby(by = 'Monomer'):
            df['Avg dist smoothed'], df['Smooth Avg dist std'] = time_smoothing(df['Avg dist'].to_list(), df['Time (ns)'].to_list(), smooth_window)
            plt.plot(df['Time (ns)'], df['Avg dist smoothed'], label = f'{mon} {vars[i]}', color=colors[color_count])
            plt.fill_between(df['Time (ns)'], df['Avg dist smoothed']-df['Smooth Avg dist std'], df['Avg dist smoothed']+df['Smooth Avg dist std'], color=colors[color_count], alpha=0.5)
            color_count += 1
    plt.xlabel('Time (ns)')
    plt.ylabel('Avg. TM distance, (' + r'$\AA$' + ')')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(0, 1000)
    #plt.ylim(0, 2)
    plt.xticks(range(0, 1001, 200))
    #plt.yticks(np.linspace(0, 2, 9))
    plt.tight_layout()
    plt.savefig(outplot_path + 'avgTMdistance_vs_time.png', dpi=600)
    plt.close()'''

    adf = pd.DataFrame.from_dict(tot_avgs)
    color_count = 0
    for v in vars:
        df = adf[adf['Variant'] == v]
        ix = df.index
        df.loc[ix, 'Avg dist smoothed'], df.loc[ix, 'Smooth Avg dist std'] = time_smoothing(df['Avg dist'].to_list(), df['Time (ns)'].to_list(), smooth_window)
        plt.plot(df['Time (ns)'], df['Avg dist smoothed'], label = f'{v}', color=colors[color_count])
        plt.fill_between(df['Time (ns)'], df['Avg dist smoothed']-df['Smooth Avg dist std'], df['Avg dist smoothed']+df['Smooth Avg dist std'], color=colors[color_count], alpha=0.5)
        color_count += 1
    plt.xlabel('Time (ns)')
    plt.ylabel('Avg. TM distance, (' + r'$\AA$' + ')')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 5)
    leg.get_frame().set_linewidth(0.0)  
    plt.xlim(0, 1000)
    #plt.ylim(0, 2)
    plt.xticks(range(0, 1001, 200))
    #plt.yticks(np.linspace(0, 2, 9))
    plt.tight_layout()
    plt.savefig(outplot_path + 'Pervariant_avgTMdistance_vs_time.png', dpi=600)
    plt.close()

def get_max_POPC_lining_freq(hole_features, monomers, out):

    results = {'Monomer':[], 'Resid':[], 'frequency%':[]}
    for i in range(len(monomers)):
        features = pd.read_csv(hole_features[i], header=0, index_col=False)
        #print(features.head())
        popcs = features[features['Resname'] == 'POP']
        #print(popcs.head(), popcs.shape)
        ntimes = len(np.unique(features['Time (ns)'].to_list()))
        for resix, df in popcs.groupby(by = 'Resid'):
            count = 0
            for t, tdf in df.groupby(by = 'Time (ns)'):
                count += 1
            results['Resid'].append(resix)
            results['frequency%'].append((100*count)/ntimes)
            results['Monomer'].append(monomers[i])


    resdf = pd.DataFrame.from_dict(results)
    resdf.to_csv(out, header = True, index=False)

vars = ['wt', 'M654V']

monomers = ['A', 'B']

for v in vars:
    struc = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/trj_data/{v}/xtcs/TMC1_{v}_p_l_noIONS.pdb'
    trj = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/trj_data/{v}/xtcs/TMC1_{v}_p_l_noIONS.xtc'
    out = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/xy-rgyr/{v}_xyrgyr.csv'
    #get_xy_gyration(struc, trj, out)

    resid = 654
    out = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/{v}_surroundings_noexclusion.csv'
    #get_surroundings(struc, trj, resid, out, exclude=0)
    out = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/Average_pos_resid{resid}_{v}.csv'
    #get_aver_position(struc, trj, resid, out, backbone=False)
    hole_features = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/hole/csv/{v}/features_{mon}_alltimes.csv' for mon in ['A', 'B']]
    final = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/{v}_chokepoints_monomer_{mon}_noexclusion.txt' for mon in ['A', 'B']]
    #compare_w_minrad(out, hole_features, final)
    out = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/POPC_freqs_{v}.csv'
    get_max_POPC_lining_freq(hole_features, monomers, out)

    out_dists = f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/othertm/{v}_S9_othertm_cog_distances.csv'
    #get_helix_z_distances(struc, trj, resid, out_dists)
    #get_helix_cog_distances(struc, trj, resid, out_dists, backbone = False, tm_helices=lastarticle_tm_helices)
outs = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/{v}_surroundings_noexclusion.csv' for v in vars]
outplot_path = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/plots/' #othertm/
#interactors_histograms(outs, vars, outplot_path)
outs = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/{v}_S9_othertm_cog_distances.csv' for v in vars]
#compare_dists_byhelix(outs, outplot_path)

outs = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/xy-rgyr/{v}_xyrgyr.csv' for v in vars]
outplot_path = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/xy-rgyr/' #othertm/
#compare_rgyrs(outs, vars, outplot_path)

outs = [f'/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/{v}_S9_othertm_cog_distances.csv' for v in vars]
outplot_path = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/plots/' #othertm/
#get_avgdists(outs, vars, outplot_path)

data = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/surroundings/no_backbone/plots/wt-M654V_interactors_diffs.csv'
#assign_regions(data)
