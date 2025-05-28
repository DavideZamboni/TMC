import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

vars = ['wt', 'M654V']
path = str(sys.argv[1])
evaluate = ['POPCs_xy']

data = {}
for var in vars:
    for mol in evaluate:
        data[f'g(r) {var} {mol}'] = []
        data['r (' + r'$\AA$' + f') {var} {mol}'] = []
        p = f'{path}/{var}/'
        with open(f'{p}rdf_TMC1_{var}_{mol}.xvg', 'r') as f:
            rel = [l for l in f.readlines() if l[0] not in ['@', '#']]
            for l in rel:
                data[f'r (' + r'$\AA$' + f') {var} {mol}'].append(float(l.split()[0])*10)
                data[f'g(r) {var} {mol}'].append(float(l.split()[1]))

for var in vars:
    vd = {}
    for k in data.keys():
        if var in k:
            vd[k] = data[k]
    df = pd.DataFrame.from_dict(vd)
    df.to_csv(f'{path}Complete_rdf_data_{var}.csv', header=True, index=False)

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
labels = {'CLA':r'$Cl^{-}$', 'POT':r'$K^{+}$', 'POPCs_noxy':'POPC', 'POPCs_xy':'POPC'}

for mol in evaluate:
    i = 0

    for var in vars:
        lab = f'{labels[mol]} {var}'
        plt.plot(data['r (' + r'$\AA$' + f') {var} {mol}'], data[f'g(r) {var} {mol}'], label = lab, color = colors[i])
        i += 1
    
    plt.xlabel('r (' + r'$\AA$)')
    plt.ylabel('Raidial Distribution Function')
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 1)
    leg.get_frame().set_linewidth(0.0)  
    plt.tight_layout()
    plt.savefig(path + f'rdf_compared_{mol}.png', dpi=600)
    plt.close()

for mol in evaluate:
    i = 0

    for var in vars:
        lab = f'{labels[mol]} {var}'
        if 'POPCs' in mol:
            rads = []
            vals = []
            for count in range(len(data['r (' + r'$\AA$' + f') {var} {mol}'])):
                if data['r (' + r'$\AA$' + f') {var} {mol}'][count] >= 10:
                    rads.append(data['r (' + r'$\AA$' + f') {var} {mol}'][count])
                    vals.append(data[f'g(r) {var} {mol}'][count])
            data['r (' + r'$\AA$' + f') {var} {mol}'] = rads
            data[f'g(r) {var} {mol}'] = vals
        plt.plot(data['r (' + r'$\AA$' + f') {var} {mol}'], data[f'g(r) {var} {mol}'], label = lab, color = colors[i])
        i += 1
    
    plt.xlabel('r (' + r'$\AA$)')
    plt.ylabel('Raidial Distribution Function')
    if 'POPCs' in mol:
        plt.xlim((10, 40))
        plt.xticks(np.linspace(10, 40, 6))
    else:
        plt.xlim((0, 40))
        plt.xticks(np.linspace(0, 40, 9))
    leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, markerscale = 1)
    leg.get_frame().set_linewidth(0.0)  
    plt.tight_layout()
    plt.savefig(path + f'rdf_compared_{mol}_reduced.png', dpi=600)
    plt.close()
