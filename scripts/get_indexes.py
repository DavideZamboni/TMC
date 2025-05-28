## To obtain GROMACS index files required for the computation of the density maps
import MDAnalysis as mda
import sys
import os
struc = str(sys.argv[1])
outpath = str(sys.argv[3])
    
reslen = 760
Tm_dict = {'S1': 'resid 180:220',
           'S2': 'resid 280:320',
           'S3': 'resid 350:380',
           'S4': 'resid 402:427',
           'S5': 'resid 435:462', 
           'S6': 'resid 517:552',
           'S7': 'resid 574:589',
           'S8': 'resid 595:620',
           'S9': 'resid 635:665',
           'S10': 'resid 700:730'}
poreforming = ['S4', 'S5', 'S6', 'S7']


## for densmaps of lipids, protein, ions

complete_idx = f'{outpath}index.ndx'
new_idx =f'{outpath}channel_index.ndx'

sele = {}
no_channel = []
channel = []
no_channel_A = []
channel_A = []
no_channel_B = []
channel_B = []
for k in Tm_dict.keys():
    resids = [int(i) for i in Tm_dict[k].lstrip('resid ').split(':')]
    if k in poreforming:
        sele[k] = f'(resid {resids[0]} to {resids[1]}) or (resid {resids[0] + reslen} to {resids[1] + reslen})'
        channel.append(f'(resid {resids[0]} to {resids[1]}) or (resid {resids[0] + reslen} to {resids[1] + reslen})')
        channel_A.append(f'(resid {resids[0]} to {resids[1]})')
        channel_B.append(f'(resid {resids[0] + reslen} to {resids[1] + reslen})')
    else:
        no_channel.append(f'(resid {resids[0]} to {resids[1]}) or (resid {resids[0] + reslen} to {resids[1] + reslen})')
        no_channel_A.append(f'(resid {resids[0]} to {resids[1]})')
        no_channel_B.append(f'(resid {resids[0] + reslen} to {resids[1] + reslen})')
sele['no_channel'] = ' or '.join(no_channel)
sele['channel'] = ' or '.join(channel)
sele['channel_CA'] = f"protein and name CA and ({sele['channel']})"
sele['POT'] = 'resname POT'
sele['CLA'] = 'resname CLA'
sele['TM'] = sele['no_channel'] + ' or ' + sele['channel']
zones = [j for j in sele.keys()]

u = mda.Universe(struc)
files = []
for k in sele.keys():
    ats = u.select_atoms(sele[k])
    fname = outpath + f'{k}.ndx'
    ats.write(fname)
    files.append(fname)

            
old_idx_lines = {}
with open(complete_idx, 'r') as f:
    ls = [l for l in f.readlines()]
    lastkey = None
    for l in ls:
        if l.startswith('[ '):
            lastkey = l.lstrip('[ ').split(' ')[0]
            old_idx_lines[lastkey] = ''
            #print(f'Reading complete_idx, found key: {lastkey}')
        elif lastkey is not None:
            old_idx_lines[lastkey] += l

with open(new_idx, 'w') as f:
    for k in old_idx_lines.keys():
        if k not in zones:
            f.write(f'[ {k} ]\n')
            f.write(old_idx_lines[k])

    for z in zones:
        with open(outpath + f'{z}.ndx', 'r') as tmp:
            mylines = [i for i in tmp.readlines() if not i.startswith('[')]
            f.write(f'[ {z} ]\n')
            #print(f'Adding zone dta from {z}.ndx.\nwill consider lines:\n{mylines}')
            f.write(''.join(mylines))

with open(new_idx, 'a') as f:
    for z in zones:
        with open(outpath + f'{z}.ndx', 'r') as tmp:
            mylines = [i for i in tmp.readlines() if not i.startswith('[')]
            #print(f'writing zone file {z}.ndx.\nwill write lines:\n{mylines}')
            f.write(f'[ {z} ]\n')
            f.write(''.join(mylines))

for z in zones:
    os.remove(outpath + f'{z}.ndx')
