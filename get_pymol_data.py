import pandas as pd
import numpy as np
from os import mkdir, remove, listdir
import MDAnalysis as mda
from json import load
from MDAnalysis.analysis import align

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def center_vector_on_onepoint(monomer_com, vector, scaling = None):
    #print(monomer_com, vector)
    if scaling is not None:
        unitary = unit_vector(vector)
        #print('unitary: ', unitary)
        #print('norm: ', np.linalg.norm(vector))
        myvector = unitary * (np.linalg.norm(vector)*scaling)
        #print('myvector: ', myvector)
        
    else:
        myvector = vector

    # Center the vector on the com
    # Compute the coordinates of the start and end points
    start_point = tuple(-v / 2 for v in myvector)
    end_point = tuple(v / 2 for v in myvector)
    #print('start_point: ', start_point)
    #print('end_point: ', end_point)

    # Calculate displacement between centres of mass, Center the vector on the point
    centered_start = [] #tuple(start_point[i] + monomer_com[i] for i in range(len(start_point)))
    centered_end = [] #tuple(end_point[i] + monomer_com[i] for i in range(len(end_point)))
    for i in range(len([j for j in monomer_com])):
        coordstart = start_point[i] + monomer_com[i]
        coordend = end_point[i] + monomer_com[i]
        centered_start.append(coordstart)
        centered_end.append(coordend)
    #print('centered_start: ', centered_start)
    #print('centered_end: ', centered_end)

    return centered_start, centered_end

def get_centroid(universe, sel, path, outsel = 'protein or resname POPC', outname = 'centroid.pdb', start = 0, stop = -1, aligned = False):
    from MDAnalysis.analysis import diffusionmap, align

    frames = [ts.frame for ts in universe.trajectory[start:stop]]
    print(frames)

    traj_folder=path+'aligned_traj/'
    try:
        mkdir(traj_folder)
    except OSError as error:
        True
    
    tstart = universe.trajectory[start].time
    tstop = universe.trajectory[stop].time
    print(f'Doing interval from frame {start} to {stop}, corresponding to times {tstart/1000} to {tstop/1000} ns')

    if not aligned:
        universe.trajectory[start]
        fr0 = f'{traj_folder}structure_frame{start}.gro'
        universe.atoms.write(fr0)
        aligned=f'{traj_folder}trajectory_alignedto_frame{start}.xtc'
        aligner = align.AlignTraj(universe, universe, select=sel, filename=aligned).run(start, stop, 1, verbose=True)
        my_uni = mda.Universe(fr0, aligned)
    else:
        print('Trajectory set as already aligned')
        my_uni = universe

    mat = diffusionmap.DistanceMatrix(my_uni, select=sel).run(0, -1, 1, verbose=True)
    matrix = mat.results.dist_matrix
    tot_rmsds = np.sum(matrix, axis = 0)
    mymin=np.argmin(tot_rmsds)
    print('Minimum RMSD AFTER alignment: ', tot_rmsds[mymin], '\nFound in frame: ', frames[mymin], ' corresponding to the structur at time ', universe.trajectory[frames[mymin]].time/1000, ' ns')

    universe.trajectory[frames[mymin]]
    out = universe.select_atoms(outsel)
    out.write(path + outname, bonds="conect")

    return universe.trajectory[frames[mymin]].time/1000

def get_dipole_data(dip_path, times, outfile, sels = ['A', 'B', 'protein']):
    outdata = {'selection':[], 'com':[], 'dip_vec':[]}
    for var in times.keys():
        for s in sels:
            dip_df = pd.read_csv(f'{dip_path}Dipoles_TMC1_{var}_{s}.csv')
            mydata = dip_df[dip_df['Time (ns)'] == times[var]]
            print(mydata.head())

def get_best_interval_selection(Tm_dict, mon_residss, monomer_number, reslen = 760):
    mon_resids = []
    for i in mon_residss:
        if i not in mon_resids:
            mon_resids.append(i)
    mon_resids.sort()
    #print(mon_resids)
    mon_intervals = []
    ordered_info = []
    #print(ordered_info)
    for k in Tm_dict.keys():
        ordered_info.append((k, Tm_dict[k][0] + reslen * monomer_number, Tm_dict[k][1] + reslen * monomer_number))
        for j in Tm_dict[k]:
            mon_intervals.append(j + reslen * monomer_number)

    ordered_info.sort(key=lambda x: int(x[1]))

    intervals = []
    for k in range(len(mon_intervals)):
        #print(k, mon_intervals[k])
        if k == 0:
            allocate = [i for i in mon_resids if i <= mon_intervals[k]]
            #print(allocate)
            if len(allocate) > 0:
                first = min(allocate)
                intervals.append(first)
            else:
                intervals.append(mon_intervals[k])

        elif k < len(mon_intervals) - 1: # not the end of the last helix
            if k % 2 == 0:  # k identifies start of helix, look which residues are between end of previous helix and this one
                allocate = [i for i in mon_resids if mon_intervals[k - 1] <= i <= mon_intervals[k]]
                if len(allocate) > 0:
                    #print(allocate)
                    diffs = [[i -mon_intervals[k - 1], mon_intervals[k] - i ] for i in allocate]
                    which_val = [i.index(min(i)) for i in diffs]
                    if 0 in which_val:
                        def_end = max([allocate[i] for i in range(len(allocate)) if which_val[i] == 0])
                        intervals.append(def_end)
                    else:
                        intervals.append(mon_intervals[k - 1])
                    
                    if 1 in which_val:
                        def_start = min([allocate[i] for i in range(len(allocate)) if which_val[i] == 1])
                        intervals.append(def_start)
                    else:
                        intervals.append(mon_intervals[k])
                else:
                    intervals.append(mon_intervals[k - 1])
                    intervals.append(mon_intervals[k])
        elif k == len(mon_intervals) - 1: # after the end of last helix
            allocate = [i for i in mon_resids if i > mon_intervals[k]]
            #print(allocate)
            if len(allocate) > 0:
                last = max(allocate)
                intervals.append(last)
            else:
                intervals.append(mon_intervals[k])

    #print(ordered_info)
    #print(mon_resids)
    #print(intervals)
    #print()
    i = 0
    extended_dict = {}
    for k in range(0, len(intervals), 2): 
        mykey = ordered_info[i][0]
        extended_dict[mykey] = (intervals[k], intervals[k+1])
        i += 1               
    return extended_dict

def get_porelining_selection(var, time, hole_csvs_folder, centroid, outdir, scripts, monomers = ['A', 'B'], around_membrane = (70, 130)):
    '''
    Writes the .pml script to get the representation o the pore-lining atoms and residues
    '''
    Tm_dict = {'S1': (176,217), 'S2': (283,303), 'S3': (354,378), 'S4': (403,428), 'S5': (436,463), 'S6': (518,553), 'S7': (575,590), 'S8': (596,618), 'S9': (637,653), 'S10': (702,728)}
    reslen = 760
    header = f'''
set depth_cue, off
set hash_max, 400
set antialias, 2
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
bg_color white
set transparency_mode, 3
set ray_shadow, 0
set ray_trace_mode, 1
set ray_trace_gain, 0
###
### PARAMS REQUIRED TO RAY-TRACE FRAMES:
set ray_trace_frames=1
set cache_frames=0
mclear
###

load {centroid}, centroid
hide everything, all

'''

    calls_df = pd.read_csv(f'{hole_csvs_folder}/spotted_intertypes.csv', header = 0, index_col = False)
    mycalls = [(atom_name, resname, call) for (atom_name, resname, call) in calls_df.itertuples(index=False)]
    identif = [(i[0], i[1]) for i in mycalls]
    #print(mycalls)
    seldict = {'Acidic': set(), 'Basic': set(), 'Hydrophobic': set(), 'Polar': set()}
    selresids = {'Acidic': set(), 'Basic': set(), 'Hydrophobic': set(), 'Polar': set()}
    popcs = set()
    backboneatoms = {'Acidic': set(), 'Basic': set(), 'Hydrophobic': set(), 'Polar': set()}
    bbresids = set()
    prot_resids = []
    backbone_names = ['N', 'C', 'CA', 'O']
    for monomer in monomers:
        datafile = f'{hole_csvs_folder}/features_{monomer}_alltimes.csv'
        res_df = pd.read_csv(datafile, header=0, usecols=['Time (ns)', 'z-coordinate (' + r'$\AA$)', 'Atom Name', 'Resname', 'Resid'])
        #print(res_df.shape)
        res_df = res_df[res_df['Time (ns)'] == time]
        #print(res_df.shape)
        if len(around_membrane) == 2:
            #print(res_df.shape)
            res_df = res_df[(res_df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (res_df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
            #print(res_df.shape)
        vals = res_df[['Atom Name', 'Resname', 'Resid']]
        dat = vals.drop_duplicates()
        
        for (atom_name, resname, resid) in dat.itertuples(index=False):
            #print(atom_name, resname, resid)
        
            if resname == 'POP':
                resname += 'C'
            else:
                prot_resids.append(resid)
            myindex = identif.index((atom_name, resname))
            call = mycalls[myindex][2]
            if atom_name not in backbone_names:
                seldict[call].add(f'(resid {resid} and resname {resname} and name {atom_name})')
                if resname != 'POPC':
                    selresids[call].add(f'(resid {resid} and resname {resname})')
                else:
                    popcs.add(f'(resid {resid} and resname {resname})')

            else:
                backboneatoms[call].add(f'(resid {resid} and resname {resname} and name {atom_name})')
                bbresids.add(f'(resid {resid} and resname {resname})')
                #print(f'{call}, (resid {resid} and resname {resname} and name {atom_name})')
    mon_resids = [[i for i in prot_resids if j*reslen <= i <= (j + 1)*reslen ] for j in range(len(monomers))]
    monomers_selections = {}
    for monomer_number in range(len(monomers)):
        mysels = get_best_interval_selection(Tm_dict, mon_resids[monomer_number], monomer_number, reslen = 760)
        selstring = ''
        for k in mysels.keys():
            selstring += f'i. {mysels[k][0]}:{mysels[k][1]} + '
        selstring = selstring[:-3]
        monomers_selections[monomers[monomer_number]] = selstring

    header += f'''
select {monomers_selections['A']}
set_name sele, TMA

select {monomers_selections['B']}
set_name sele, TMB

set cartoon_cylindrical_helices, 1

hide everything, resn POPC

show cartoon, TMA
show cartoon, TMB
set cartoon_transparency, 0.4

color black, TMA
color gray70, TMB
rebuild

'''
    for k in seldict.keys():
        header += f'select {" or ".join(seldict[k])}\nset_name sele, {k}atoms\n'
    
    for k in selresids.keys():
        header += f'select {" or ".join(selresids[k])}\nset_name sele, {k}res\n'

    header += f'select {" or ".join(popcs)}\nset_name sele, POPCres\n'

    header += f'select {" or ".join(bbresids)}\ncreate mybbs, sele\nhide everything, mybbs\n'

    for k in backboneatoms.keys():
        if len(backboneatoms[k]) > 0:
            header += f'select mybbs and {" or ".join(backboneatoms[k])}\nset_name sele, {k}_bb_atoms\n'


    header += f'''
as sticks, (Acidicres or Polarres or Hydrophobicres or POPCres or mybbs) and not (hydrogen or backbone)
as sticks, Basicres and not (hydrogen or backbone)

color paleyellow, POPCres
color olive, Hydrophobicatoms
color olive, Hydrophobic_bb_atoms
color skyblue, Polaratoms
color skyblue, Polar_bb_atoms
color hotpink, Basicatoms
color hotpink, Basic_bb_atoms
color limegreen, Acidicatoms
color limegreen, Acidic_bb_atoms

set stick_radius, 0.6
set stick_radius, 0.2, POPCres
as spheres, Acidicatoms or Polaratoms or Hydrophobicatoms
as spheres, Basicatoms
as spheres, Acidic_bb_atoms or Polar_bb_atoms or Hydrophobic_bb_atoms
as spheres, Basic_bb_atoms

load {scripts}cgo_grid.py
cgo_grid pos1=[0, 0, 120], pos2=[0, 150, 120], pos3=[150, 0, 120], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Extracellular"
cgo_grid pos1=[0, 0, 80], pos2=[0, 150, 85], pos3=[150, 0, 80], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Cytosolic"

set cgo_transparency, 0.5, Extracellular
set cgo_transparency, 0.5, Cytosolic
set sphere_quality, 3
bg_color white


viewport 800, 600
refresh

mset 1 x741
## MEMBRANE

set_view (\
    -0.999684513,   -0.018730300,   -0.016624127,\
    -0.016235318,   -0.020742301,    0.999652267,\
    -0.019068565,    0.999609172,    0.020431632,\
     0.000173620,   -0.000122726, -392.634033203,\
    71.920883179,   69.163581848,   96.978775024,\
   266.117340088,  519.149047852,  -20.000000000 )
  
refresh
mview store, 1
mview reinterpolate
refresh

turn y, 120
mview store, 61, power=1.0
mview reinterpolate
turn y, 120
mview store, 121, power=1.0
mview reinterpolate
turn y, 120
mview store, 181, power=1.0
mview reinterpolate
refresh

## CHAN A

set_view (\
    -0.999684513,   -0.018730300,   -0.016624127,\
    -0.016235318,   -0.020742301,    0.999652267,\
    -0.019068565,    0.999609172,    0.020431632,\
     0.000089571,   -0.000053971, -221.500350952,\
    27.978700638,   68.394294739,   97.731880188,\
    200,  250,  -20.000000000 )
refresh

mview store, 231
mview reinterpolate
refresh

turn y, 120
mview store, 291, power=1.0
mview reinterpolate
turn y, 120
mview store, 351, power=1.0
mview reinterpolate
turn y, 120
mview store, 411, power=1.0
mview reinterpolate
refresh

## MEMBRANE

set_view (\
    -0.999684513,   -0.018730300,   -0.016624127,\
    -0.016235318,   -0.020742301,    0.999652267,\
    -0.019068565,    0.999609172,    0.020431632,\
     0.000173620,   -0.000122726, -392.634033203,\
    71.920883179,   69.163581848,   96.978775024,\
   266.117340088,  519.149047852,  -20.000000000 )
  
refresh
mview store, 461
mview reinterpolate
refresh

## CHAN B

set_view (\
    -0.997482836,    0.023012472,    0.067043528,\
     0.065815561,   -0.050500486,    0.996552050,\
     0.026318962,    0.998458385,    0.048858926,\
     0.000323340,   -0.000051846, -219.115539551,\
    92.045722961,   69.517112732,   95.743751526,\
    200,  250,  -20.000000000 )
    
refresh

mview store, 511
mview reinterpolate
refresh

turn y, 120
mview store, 571, power=1.0
mview reinterpolate
turn y, 120
mview store, 631, power=1.0
mview reinterpolate
turn y, 120
mview store, 691, power=1.0
mview reinterpolate
refresh

refresh

## MEMBRANE

set_view (\
    -0.999684513,   -0.018730300,   -0.016624127,\
    -0.016235318,   -0.020742301,    0.999652267,\
    -0.019068565,    0.999609172,    0.020431632,\
     0.000173620,   -0.000122726, -392.634033203,\
    71.920883179,   69.163581848,   96.978775024,\
   266.117340088,  519.149047852,  -20.000000000 )
  
refresh
mview store, 741
mview reinterpolate
refresh
#mpng {outdir}/movie_pore_lining/{var}/{var}_f, 1, 741, width=3300, height=2400 # Uncomment to perorm rendering of the movie
    '''
    with open(f'{outdir}/movie_pore_lining_{var}.pml', 'w') as f:
        f.write(header)

def get_VdW_pml(var, centroid, outdir, scripts):
    header = f'''
### QuteMol-like params
set depth_cue, off
set hash_max, 24000
set antialias, 2
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
bg_color white
set transparency_mode, 3
set ray_shadow, 0
###

### PARAMS REQUIRED TO RAY-TRACE FRAMES:
set ray_trace_frames=1
set cache_frames=0
mclear
###

load {outdir}HOLE_centers_centroid_{var}.pdb, cens
load {centroid}, centroid

select chain A
set_name sele, monA

select chain B
set_name sele, monB

select resid 403-428 or resid 1163-1188
set_name sele, S4
select resid 436-463 or resid 1196-1223
set_name sele, S5
select resid 518-553 or resid 1278-1313
set_name sele, S6
select resid 575-590 or resid 1335-1350
set_name sele, S7
select resid 1 or resid 761
set_name sele, Nter
select resid 760 or resid 1520
set_name sele, Cter

set cartoon_cylindrical_helices, 1
hide everything, resn POPC
show cartoon, monA
show cartoon, monB
set cartoon_transparency, 0.4
show spheres, Nter
show spheres, Cter
alter Nter, vdw=0.5
alter Cter, vdw=0.5
rebuild

color black, monA
color gray70, monB
color blue, S4
color green, S5
color purple, S6
color red, S7

load {scripts}cgo_grid.py
cgo_grid pos1=[0, 0, 120], pos2=[0, 150, 120], pos3=[150, 0, 120], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Extracellular"
cgo_grid pos1=[0, 0, 85], pos2=[0, 150, 85], pos3=[150, 0, 85], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Cytosolic"

select chain C
set_name sele, centers
create centers2, centers
alter centers2, vdw=0.3
alter centers, vdw=b
rebuild

show mesh, centers
color tv_orange, centers
color hotpink, centers2
hide spheres, centers
refresh

viewport 800, 600
refresh

mset 1 x741
## MEMBRANE

set_view (\
    -0.989879668,   -0.002948569,    0.141867757,\
     0.141889721,   -0.009443900,    0.989835501,\
    -0.001578606,    0.999951005,    0.009766647,\
     0.000173620,   -0.000122726, -363.140228271,\
    71.920883179,   69.163581848,   96.978775024,\
   236.623641968,  489.655364990,  -20.000000000 )

refresh
mview store, 1
mview reinterpolate
refresh

turn y, 120
mview store, 61, power=1.0
mview reinterpolate
turn y, 120
mview store, 121, power=1.0
mview reinterpolate
turn y, 120
mview store, 181, power=1.0
mview reinterpolate
refresh

## CHAN A

set_view (\
    -0.999684513,   -0.018730300,   -0.016624127,\
    -0.016235318,   -0.020742301,    0.999652267,\
    -0.019068565,    0.999609172,    0.020431632,\
     0.000089571,   -0.000053971, -221.500350952,\
    27.978700638,   68.394294739,   97.731880188,\
    200,  250,  -20.000000000 )
refresh

mview store, 231
mview reinterpolate
refresh

turn y, 120
mview store, 291, power=1.0
mview reinterpolate
turn y, 120
mview store, 351, power=1.0
mview reinterpolate
turn y, 120
mview store, 411, power=1.0
mview reinterpolate
refresh

## MEMBRANE

set_view (\
    -0.989879668,   -0.002948569,    0.141867757,\
     0.141889721,   -0.009443900,    0.989835501,\
    -0.001578606,    0.999951005,    0.009766647,\
     0.000173620,   -0.000122726, -363.140228271,\
    71.920883179,   69.163581848,   96.978775024,\
   236.623641968,  489.655364990,  -20.000000000 )
  
refresh
mview store, 461
mview reinterpolate
refresh

## CHAN B

set_view (\
    -0.997482836,    0.023012472,    0.067043528,\
     0.065815561,   -0.050500486,    0.996552050,\
     0.026318962,    0.998458385,    0.048858926,\
     0.000323340,   -0.000051846, -219.115539551,\
    92.045722961,   69.517112732,   95.743751526,\
    200,  250,  -20.000000000 )
    
refresh

mview store, 511
mview reinterpolate
refresh

turn y, 120
mview store, 571, power=1.0
mview reinterpolate
turn y, 120
mview store, 631, power=1.0
mview reinterpolate
turn y, 120
mview store, 691, power=1.0
mview reinterpolate
refresh

refresh

## MEMBRANE

set_view (\
    -0.989879668,   -0.002948569,    0.141867757,\
     0.141889721,   -0.009443900,    0.989835501,\
    -0.001578606,    0.999951005,    0.009766647,\
     0.000173620,   -0.000122726, -363.140228271,\
    71.920883179,   69.163581848,   96.978775024,\
   236.623641968,  489.655364990,  -20.000000000 )

refresh
mview store, 741
mview reinterpolate
refresh
#mpng {outdir}/movie_vdw/{var}/{var}_f, 1, 741, width=3300, height=2400

'''
    with open(outdir + f'movie_vdw_{var}.pml', 'w') as f:
        f.write(header)

def get_VdW_channels(var, time, hole_csvs_folder, centroid, outdir, scripts, monomers = ['A', 'B'], around_membrane = (70, 130)):
    centers = ''
    myid = mda.Universe(centroid).atoms.ids[-1] + 1
    myresid = mda.Universe(centroid).atoms.resids[-1] + 1

    for monomer in monomers:
        datafile = f'{hole_csvs_folder}/features_{monomer}_alltimes.csv'
        res_df = pd.read_csv(datafile, header=0, usecols=['Time (ns)', 'x-coordinate (' + r'$\AA$)', 'y-coordinate (' + r'$\AA$)' , 'z-coordinate (' + r'$\AA$)', 'Radius ' + r'$\AA$' ])
        #print(res_df.shape)
        res_df = res_df[res_df['Time (ns)'] == time]
        #print(res_df.shape)
        if len(around_membrane) == 2:
            #print(res_df.shape)
            res_df = res_df[(res_df['z-coordinate (' + r'$\AA$)'] >= around_membrane[0]) & (res_df['z-coordinate (' + r'$\AA$)'] <= around_membrane[1])]
            #print(res_df.shape)
        res_df = res_df[[ 'x-coordinate (' + r'$\AA$)', 'y-coordinate (' + r'$\AA$)' , 'z-coordinate (' + r'$\AA$)', 'Radius ' + r'$\AA$' ]]
        for (x, y, z, r) in res_df.itertuples(index = False):
            centers += f'ATOM  {myid:>5}  CA  CEN C{myresid:>4}     {x:>7.3f} {y:>7.3f} {z:>7.3f}  1.00 {r:>5.2f}      C    X\n'
            myid += 1
            myresid += 1
            if myid >= 100000:
                myid = 1
            if myresid >= 10000:
                myresid = 1
    
    with open(outdir + f'HOLE_centers_centroid_{var}.pdb', 'w') as f:
        f.write(centers)

    get_VdW_pml(var, centroid, outdir, scripts)

def get_dipoles(var, time, dipoles_csvs_folder, centroid, outdir, scripts, scaling = 0.3, rgbs =  ("0.3, 0.3, 1.0"	, "1.0, 1.0, 0.2" , "0.2, 1.0, 0.2"), colors=['tv_blue', 'tv_yellow', 'tv_green'], selections = ['chainID A', 'chainID B', 'protein']):
    header = f'''

### QuteMol-like params
set depth_cue, off
set hash_max, 24000
set antialias, 2
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
bg_color white
set transparency_mode, 3
set ray_shadow, 0
###
### PARAMS REQUIRED TO RAY-TRACE FRAMES:
set ray_trace_frames=1
set cache_frames=0
mclear
###

load {centroid}, centroid

select chain A
set_name sele, monA

select chain B

set_name sele, monB

select resid 1 or resid 761
set_name sele, Nter

select resid 760 or resid 1520
set_name sele, Cter

hide everything, resn POPC

show cartoon, monA
show cartoon, monB
set cartoon_transparency, 0.4

select resn POPC & (e. O+P+N)
set_name sele, POPCHs
show spheres, POPCHs
show spheres, Nter
show spheres, Cter

alter Nter, vdw=0.5
alter Cter, vdw=0.5
alter POPCHs, vdw=0.5

rebuild

color tv_orange, monA
color lightmagenta, monB
color black, POPCHs
set sphere_transparency, 0.4, POPCHs

python
from __future__ import print_function
from pymol.cgo import *    # get constants
from math import *
from pymol import cmd

def modevectors(start_point, end_point, outname="vector", head=1.0, tail=0.3, head_length=1.5, headrgb="1.0,1.0,1.0", tailrgb="1.0,1.0,1.0", cutoff=4.0, skip=0, cut=0.5, atom="CA", stat="show", factor=1.0, notail=0):

    objectname = outname
    factor = float(factor)
    arrow_head_radius = float(head)
    arrow_tail_radius = float(tail)
    arrow_head_length = float(head_length)
    cutoff = float(cutoff)
    skip = int(skip)
    cut = float(cut)
    atomtype = atom.strip('"[]()')
    objectname = objectname.strip('"[]()')

    headrgb = headrgb.strip('" []()')
    tailrgb = tailrgb.strip('" []()')
    print(headrgb, tailrgb)
    hr, hg, hb = list(map(float, headrgb.split(',')))
    tr, tg, tb = list(map(float, tailrgb.split(',')))
    print(hr, hg, hb, tr, tg, tb)

    version = cmd.get_version()[1]
    arrow = []
    arrowhead = []
    arrowtail = []
    x1 = start_point[0]
    y1 = start_point[1]
    z1 = start_point[2]
    x2 = end_point[0]
    y2 = end_point[1]
    z2 = end_point[2]

    save_view = cmd.get_view(output=1, quiet=1)

    arrowtail = []

    vectorx = x2 - x1
    vectory = y2 - y1
    vectorz = z2 - z1
    length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)

    t = 1.0 - (cut / length)
    length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)
    d = arrow_head_length  # Distance from arrow tip to arrow base
    t = 1.0 - (d / length)
    if notail:
        t = 0
    tail = [
        # Tail of cylinder
        CYLINDER, x1, y1, z1\
        , x1 + (t + 0.01) * vectorx, y1 + (t + 0.01) * vectory, z1 + (t + 0.01) * vectorz\
        , arrow_tail_radius, tr, tg, tb, tr, tg, tb  # Radius and RGB for each cylinder tail
    ]
    if notail == 0:
        arrow.extend(tail)
    x = x1 + t * vectorx
    y = y1 + t * vectory
    z = z1 + t * vectorz
    dx = x2 - x
    dy = y2 - y
    dz = z2 - z
    seg = d / 100
    intfactor = int(factor)
    if version < 1.1:  # Version >= 1.1 has cone primitive
        for i in range(100, 0, -1):  # i=100 is tip of cone
            #print(i)
            t1 = seg * i
            t2 = seg * (i + 1)
            radius = arrow_head_radius * (1.0 - i / (100.0))  # Radius of each disc that forms cone
            head = [
                CYLINDER, x + t2 * dx, y + t2 * dy, z + t2 * dz\
                , x + t1 * dx, y + t1 * dy, z + t1 * dz\
                , radius, hr, hg, hb, hr, hg, hb  # Radius and RGB for slice of arrow head
            ]
            arrow.extend(head)
    else:
        head = [
            CONE, x, y, z, x + d * dx, y + d * dy, z + d * dz, arrow_head_radius, 0.0, hr, hg, hb, hr, hg, hb, 1.0, 1.0]
        arrow.extend(head)

    cmd.bg_color(color="white")
    cmd.set_view(save_view)
    cmd.load_cgo(arrow, objectname)

    return
'''
    data = {}
    for s in selections:
        dipoles = pd.read_csv(f'{dipoles_csvs_folder}/Dipoles_TMC1_{var}_{"_".join(s.split())}.csv', header = 0, index_col= False)
        dipoles = dipoles[dipoles['Time (ns)'] == time]
        #print(s, dipoles.head())
        vector = (float(dipoles['Dip_x']), float(dipoles['Dip_y']), float(dipoles['Dip_z']))
        monomer_com = (float(dipoles['com_x']), float(dipoles['com_y']), float(dipoles['com_z']))
        start, end = center_vector_on_onepoint(monomer_com, vector, scaling = scaling)
        data[s] = (monomer_com, start, end)
    print(data)
    for i in range(len(selections)):
        header += f'modevectors({data[selections[i]][1]}, {data[selections[i]][2]}, outname="{"_".join(selections[i].split())}_dip", headrgb="{rgbs[i]}", tailrgb="{rgbs[i]}", head=1.5, tail=0.5)\n'
        header += f'cmd.pseudoatom("{selections[i][-1]}_com", pos={data[selections[i]][0]}, name="COM{selections[i][-1]}")\n'
    
    header += 'python end\n'

    k = 0
    for i in selections:
        header += f'alter {i[-1]}_com, vdw=1.5\nrebuild\nshow spheres, {i[-1]}_com\ncolor {colors[k]}, {i[-1]}_com\n'
        header += f'set stick_radius, 1.5, {"_".join(i.split())}_dip\n'
        k += 1

    header += f'''

set sphere_quality, 3
bg_color white

viewport 800, 600
refresh

mset 1 x381

## MEMBRANE
set_view (\
    -0.999684572,   -0.018730294,   -0.016623914,\
    -0.016235143,   -0.020742295,    0.999652326,\
    -0.019068563,    0.999608994,    0.020431640,\
     0.000173620,   -0.000122726, -470.604339600,\
    71.920883179,   69.163581848,   96.978775024,\
  -7150.730468750, 8091.936035156,  -19.999998093 )

  
mview store, 1
mview reinterpolate
refresh

turn y, 120
mview store, 61, power=1.0
mview reinterpolate
turn y, 120
mview store, 121, power=1.0
mview reinterpolate
turn y, 120
mview store, 181, power=1.0
mview reinterpolate
refresh

## EXTRACELLULAR

set_view (\
    -0.994113803,    0.040512979,   -0.100466251,\
    -0.048039418,   -0.996125400,    0.073657051,\
    -0.097093552,    0.078049965,    0.992211163,\
     0.000099355,    0.000054598, -476.033782959,\
    66.836647034,   79.688079834,   96.122123718,\
  -12353.498046875, 13305.693359375,  -20.000000000 )

mview store, 231, power=1.0
mview reinterpolate
refresh

## MEMBRANE SHORT
set_view (\
    -0.999684572,   -0.018730294,   -0.016623914,\
    -0.016235143,   -0.020742295,    0.999652326,\
    -0.019068563,    0.999608994,    0.020431640,\
     0.000173620,   -0.000122726, -470.604339600,\
    71.920883179,   69.163581848,   96.978775024,\
  -7150.730468750, 8091.936035156,  -19.999998093 )

mview store, 281, power=1.0
mview reinterpolate
refresh

## CYTOSOLIC
set_view (\
    -0.989254177,   -0.065073714,   -0.130926818,\
    -0.051745489,    0.993361771,   -0.102750495,\
     0.136744201,   -0.094871700,   -0.986053228,\
    -0.000020146,    0.000206083, -513.443176270,\
    72.554748535,   78.198410034,  137.145507812,\
   386.913970947,  639.945739746,  -20.000000000 )

mview store, 331, power=1.0
refresh

## MEMBRANE END
set_view (\
    -0.999684572,   -0.018730294,   -0.016623914,\
    -0.016235143,   -0.020742295,    0.999652326,\
    -0.019068563,    0.999608994,    0.020431640,\
     0.000173620,   -0.000122726, -470.604339600,\
    71.920883179,   69.163581848,   96.978775024,\
  -7150.730468750, 8091.936035156,  -19.999998093 )


mview store, 381, power=1.0
mview reinterpolate
refresh

#mpng {outdir}/movie_dipoles/{var}/{var}_f, 1, 381, width=3300, height=2400
'''
    with open(outdir + f'movie_dipoles_{var}.pml', 'w') as f:
        f.write(header)

def get_ions_pml(var, centroid, outdir, scripts, dist = 50):
    ions_centroid = centroid[:-4] + '_slowdiff_ions.pdb'
    mycom = mda.Universe(centroid).select_atoms('protein').center_of_mass()

    header = f'''
### QuteMol-like params
set depth_cue, off
set hash_max, 24000
set antialias, 2
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
bg_color white
set transparency_mode, 3
set ray_shadow, 0

###

### PARAMS REQUIRED TO RAY-TRACE FRAMES:
set ray_trace_frames=1
set cache_frames=0
mclear
###

load {ions_centroid}, centroid

select chain A
set_name sele, monA

select chain B
set_name sele, monB

select resid 403-428 or resid 1163-1188
set_name sele, S4

select resid 436-463 or resid 1196-1223
set_name sele, S5

select resid 518-553 or resid 1278-1313
set_name sele, S6

select resid 575-590 or resid 1335-1350
set_name sele, S7

select resid 1 or resid 761
set_name sele, Nter

select resid 760 or resid 1520
set_name sele, Cter

set cartoon_cylindrical_helices, 1

hide everything, resn POPC

show cartoon, monA
show cartoon, monB
set cartoon_transparency, 0.4

show spheres, Nter
show spheres, Cter

alter Nter, vdw=0.5
alter Cter, vdw=0.5

select (i. 81:130 + i. 841:890) & ! e. h
set_name sele, weakCIB

select (i. 305:344 + i. 1065:1104) & ! e. h
set_name sele, strongCIB

rebuild

color black, monA
color gray70, monB
color blue, S4
color green, S5
color purple, S6
color red, S7

load {scripts}cgo_grid.py
cgo_grid pos1=[0, 0, 120], pos2=[0, 150, 120], pos3=[150, 0, 120], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Extracellular"
cgo_grid pos1=[0, 0, 80], pos2=[0, 150, 80], pos3=[150, 0, 80], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Cytosolic"

set cgo_transparency, 0.5, Extracellular
set cgo_transparency, 0.5, Cytosolic

python

from pymol import cmd

coordsTot  = [{mycom[0]},{mycom[1]},{mycom[2]}]
xcT = coordsTot[0]
ycT = coordsTot[1]
zcT = coordsTot[2]
cmd.pseudoatom('com', pos=[xcT, ycT, zcT], name='CA')

python end

hide everything, resn CLA or resn POT

select (resn CLA or resn POT) within {dist} of com
set_name sele, INW

select INW and resn CLA
set_name sele, INC

select INW and resn POT
set_name sele, INP

color cyan, INW
color tv_orange, INC

show spheres, INW

alter INW, vdw=1

rebuild

set sphere_quality, 3
bg_color white

#viewport 800, 600
refresh

mset 1 x381

## MEMBRANE

set_view (\
    -0.999492764,    0.027088631,   -0.016624127,\
    -0.017168539,   -0.019976716,    0.999652267,\
     0.026747271,    0.999433160,    0.020431632,\
     0.000173620,   -0.000122726, -366.475311279,\
    71.920883179,   69.163581848,   96.978775024,\
   239.958511353,  492.990203857,  -20.000000000 )

mview store, 1
mview reinterpolate
refresh

turn y, 120
mview store, 61, power=1.0
mview reinterpolate
turn y, 120
mview store, 121, power=1.0
mview reinterpolate
turn y, 120
mview store, 181, power=1.0
mview reinterpolate
refresh

## EXTRACELLULAR

set_view (\
    -0.997444749,    0.070977829,   -0.007886893,\
    -0.071414910,   -0.991720140,    0.106725842,\
    -0.000246601,    0.107016720,    0.994258165,\
     0.000173620,   -0.000122726, -277.906616211,\
    71.920883179,   69.163581848,   96.978775024,\
   151.389587402,  404.421386719,  -20.000000000 )

mview store, 231, power=1.0
mview reinterpolate
refresh

## MEMBRANE SHORT
set_view (\
    -0.999492764,    0.027088631,   -0.016624127,\
    -0.017168539,   -0.019976716,    0.999652267,\
     0.026747271,    0.999433160,    0.020431632,\
     0.000173620,   -0.000122726, -366.475311279,\
    71.920883179,   69.163581848,   96.978775024,\
   239.958511353,  492.990203857,  -20.000000000 )

mview store, 281, power=1.0
mview reinterpolate
refresh

## CYTOSOLIC
set_view (\
    -0.975626826,   -0.204683557,   -0.079099394,\
    -0.184357360,    0.960066199,   -0.210447177,\
     0.119015940,   -0.190735519,   -0.974400282,\
    -0.000027880,    0.000187710, -333.803405762,\
    76.277526855,   82.172615051,  134.534301758,\
   207.280578613,  460.312194824,  -20.000000000 )

mclear
mview store, 331, power=1.0
mview reinterpolate
refresh


## MEMBRANE END
set_view (\
    -0.999492764,    0.027088631,   -0.016624127,\
    -0.017168539,   -0.019976716,    0.999652267,\
     0.026747271,    0.999433160,    0.020431632,\
     0.000173620,   -0.000122726, -366.475311279,\
    71.920883179,   69.163581848,   96.978775024,\
   239.958511353,  492.990203857,  -20.000000000 )


mview store, 381
mview reinterpolate
refresh

#mpng {outdir}/movie_ions/{var}/{var}_f, 1, 381, width=3300, height=2400
'''
    with open(outdir + f'movie_ions_{var}.pml', 'w') as f:
        f.write(header)

def get_slowdiff_ions(var, centroid, outdir, scripts, ions = [('chloride', 'CLA', 'CL'), ('potassium', 'POT', 'K')], filter = False, dist = 50):
    pdb_ions = ''
    myid = mda.Universe(centroid).atoms.ids[-1] + 1
    myresid = mda.Universe(centroid).atoms.resids[-1] + 1
    positions = {1: [1], 2: [1, 2], 3: [1, 2, 3], 4: [0, 1, 2, 3]}

    oname = centroid[:-4] + '_slowdiff_ions.pdb'
    for ion in ions:
        if not filter:
            pos = pd.read_csv(f'{outdir}/filtered_lowD_{ion[0]}_{var}.csv', header=0, index_col=False, usecols = ['x', 'y', 'z'])
            namelen = positions[len(ion)]
            split = [ion[1][namelen.index(i)] if i in namelen else ' ' for i in range(4)]
            myname = ''.join(split)
            for (x, y, z) in pos.itertuples(index = False):
                pdb_ions += f'ATOM  {myid:>5} {myname}{myname} I{myresid:>4}     {x:>7.3f} {y:>7.3f} {z:>7.3f}  1.00  0.00      I   {ion[2]:>2}\n'
                myid += 1
                myresid += 1
                if myid >= 100000:
                    myid = 1
                if myresid >= 10000:
                    myresid = 1

    atoms = []
    conects = []
    with open(centroid, 'r') as f:
        l = f.readlines()
        for line in l:
            if not line.startswith('END') and not line.startswith('CONECT'):
                atoms.append(line)
            else:
                conects.append(line)
    
    with open(oname, 'w') as f:
        for s in atoms:
            f.write(s)
        for s in pdb_ions:
            f.write(s)
        for s in conects:
            f.write(s)

    get_ions_pml(var, centroid, outdir, scripts, dist = dist)

def get_CIB2_pml(var, centroid, outdir, scripts):
    header = f'''
### QuteMol-like params
set depth_cue, off
set hash_max, 24000
set antialias, 2
set light_count,8
set spec_count,1
set shininess, 10
set specular, 0.25
set ambient,0
set direct,0
set reflect,1.5
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
bg_color white
set transparency_mode, 3
set ray_shadow, 0

### PARAMS REQUIRED TO RAY-TRACE FRAMES:
set ray_trace_frames=1
set cache_frames=0
mclear
###

load {centroid}, centroid

select chain A
set_name sele, monA

select chain B
set_name sele, monB

select resid 403-428 or resid 1163-1188
set_name sele, S4

select resid 436-463 or resid 1196-1223
set_name sele, S5

select resid 518-553 or resid 1278-1313
set_name sele, S6

select resid 575-590 or resid 1335-1350
set_name sele, S7

select resid 1 or resid 761
set_name sele, Nter

select resid 760 or resid 1520
set_name sele, Cter

set cartoon_cylindrical_helices, 1

hide everything, resn POPC

show cartoon, monA
show cartoon, monB
set cartoon_transparency, 0.4

show spheres, Nter
show spheres, Cter

alter Nter, vdw=0.5
alter Cter, vdw=0.5

select (i. 86:136 + i. 846:896) & ! e. h
set_name sele, weakCIB

select (i. 304:358 + i. 1064:1118) & ! e. h
set_name sele, strongCIB

rebuild

color black, monA
color gray70, monB
#color blue, S4 better without helices colorcoding 
#color green, S5 
#color purple, S6 
#color red, S7 

load {scripts}/cgo_grid.py
cgo_grid pos1=[0, 0, 120], pos2=[0, 150, 120], pos3=[150, 0, 120], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Extracellular"
cgo_grid pos1=[0, 0, 85], pos2=[0, 150, 85], pos3=[150, 0, 85], length_x=150, length_z=150, npoints_x=30, npoints_z=30, thickness=2, color=[0.75, 0.75, 0.75], mode=1, name="Cytosolic"

set cgo_transparency, 0.5, Extracellular
set cgo_transparency, 0.5, Cytosolic

rebuild

set sphere_quality, 3
bg_color white

#viewport 800, 600
refresh

## CYTOSOLIC
set_view (\
    -0.975626826,   -0.204683557,   -0.079099394,\
    -0.184357360,    0.960066199,   -0.210447177,\
     0.119015940,   -0.190735519,   -0.974400282,\
    -0.000027880,    0.000187710, -333.803405762,\
    76.277526855,   82.172615051,  134.534301758,\
   207.280578613,  460.312194824,  -20.000000000 )


## CYTOSOLIC CIB2
# here as sticks
color tv_red, weakCIB
color marine, strongCIB

# Sticks
as sticks, weakCIB
as sticks, strongCIB
hide everything, resname CLA
hide everything, resname POT
refresh

ray 3300, 2400
png {outdir}/movie_ions/{var}/CIB2_sites_{var}.png

'''    
    with open(outdir + f'CIB2_frames_{var}.pml', 'w') as f:
        f.write(header)

def get_pymol_representations(var, time, hole_csvs_folder, dipoles_csvs_folder, centroid, outdir, scripts, monomers = ['A', 'B'], around_membrane = (70, 130), ions = [('chloride', 'CLA', 'CL'), ('potassium', 'POT', 'K')], dist = 50):
    '''
    To create the various folders to store the frames in order to make the movies.
    From the output directory, it will create 4 folders:
    - movie_pore_lining
    - movie_vdw
    - movie_ions
    - movie_dipoles
    Each will have a separate folder for each variant.
    '''
    folders = ['movie_pore_lining', 'movie_vdw', 'movie_ions', 'movie_dipoles']
    for f in folders:
        try:
            mkdir(f'{outdir}/{f}/')
        except:
            True

        try:
            mkdir(f'{outdir}/{f}/{var}/')
        except:
            True

    get_slowdiff_ions(var, centroid, outdir, scripts, ions = ions)
    get_CIB2_pml(var, centroid, outdir, scripts)
    get_porelining_selection(var, time, hole_csvs_folder, centroid, outdir, scripts, monomers = monomers, around_membrane = around_membrane)
    get_VdW_channels(var, time, hole_csvs_folder, centroid, outdir, scripts, monomers = ['A', 'B'], around_membrane = (70, 130))
    get_dipoles(var, time, dipoles_csvs_folder, centroid, outdir, scripts, selections = ['chainID A', 'chainID B', 'protein'])


def calculate_instantaneous_diffusion(df):
    # Get the initial and final positions
    #print(df.iloc[:, 1:-1])
    diffs = []
    permuts = []
    for i in range(df.shape[0]):
        for j in range(i + 1, df.shape[0]):
            permuts.append((i, j))

    for (i, j) in permuts:
        initial_pos = df.iloc[i, 1:-1].values
        final_pos = df.iloc[j, 1:-1].values

        # Calculate displacement vector
        displacement = final_pos - initial_pos

        # Calculate squared displacement magnitude
        squared_displacement = np.sum(displacement ** 2)

        # Estimate diffusion coefficient using Einstein's relation
        time_difference = df.iloc[j]['Time (ns)'] - df.iloc[i]['Time (ns)']
        #print(time_difference)

        ## NW: diffusion_coefficient in nm**2/ps
        ## the thresholds I get from GROMACS msd files are in 1e-5 cm**2/s = 1e-9 m**2/s
        ## need to convert to the right units, here  in nm**2/ps
        ## 1 nm = 1e-9 m, 1 ps = 1e-12 s
        ## 1 nm**2 = 1e-18 m**2
        ## 1 nm**2 / ps = 1e-18 m**2 / 1e-12 s = 1e-6 m**2/s ==> need to multiply the instantaneous diffusion coefficient by 1e-3
        conv_factor = 1e-3
        diffusion_coefficient = (squared_displacement / (6 * time_difference)) * conv_factor

        #print(i, j, diffusion_coefficient)
        diffs.append(diffusion_coefficient)
    #print(np.mean(diffs))
    return diffs

def all_lower(diffs, threshold):
    lower = True
    for d in diffs:
        if d > threshold:
            lower = False
    return lower

def filter_positions(lowD_file, threshold, tstep = 0.1):
    df = pd.read_csv(lowD_file, header=0, index_col=False)
    data = {'x':[], 'y':[], 'z':[], 'ion':[]}
    for ion, idf in df.groupby(by = 'loc'):
        idf.sort_values(by = 'Time (ns)', ignore_index = True, inplace = True)
        times =  [t for t in idf['Time (ns)']]
        gaps = [i for  i in range(1, len(times)) if abs(times[i] - times[i-1]) > (tstep * 1.5)]
        if len(gaps) < 0:
            gaps = [idf.shape[0]]
        start = 0
        last = 0
        while last < len(gaps):
            for i in range(start+1, gaps[last], 2):
                if i + 1 >= gaps[last]:
                    points = idf.iloc[i-2:i+1]
                else:
                    points = idf.iloc[i-1:i+2]

                diffusions = calculate_instantaneous_diffusion(points)
                if all_lower(diffusions, threshold):
                    lowpos = points.iloc[1]
                    #print(points)
                    #print(lowpos['x'],lowpos['y'],lowpos['z'])
                    data['x'].append(lowpos['x'])
                    data['y'].append(lowpos['y'])
                    data['z'].append(lowpos['z'])
                    data['ion'].append(ion)
                #print(points)
                #print(start, last, i)            
            start = gaps[last]
            last += 1
    return pd.DataFrame.from_dict(data)

Tm_dict = {'S1': 'resid 176:217', 'S2': 'resid 283:303', 'S3': 'resid 354:378', 'S4': 'resid 403:428', 'S5': 'resid 436:463', 'S6': 'resid 518:553', 'S7': 'resid 575:590', 'S8': 'resid 596:618', 'S9': 'resid 637:653', 'S10': 'resid 702:728'}
channel_helices = ['S4', 'S5', 'S6', 'S7']
monomers = ['A', 'B']
vars = ['wt']
index = {'wt':12344, 'M654V':12343}
times = {}
basepath = '/home/davide/storage/lab_dati/davide_TMC1/definitivi/'
renderings = f'{basepath}/renderings/'
try:
    mkdir(renderings)
except:
    True

## Get the centroid, LONG
'''
for var in vars:
    path = f'{basepath}trj_data/{var}/xtcs/'
    universe = mda.Universe(f'{path}TMC1_{var}_p_l.pdb', f'{path}TMC1_{var}_p_l.xtc')
    tmres = '('
    for i in Tm_dict.keys():
        tmres += Tm_dict[i] + ' or '
    tmres = tmres[:-4] + ')'
    idx = index[var]
    ch = 'backbone and ('
    for mon in range(len(monomers)):
        if mon == 0:
            ch += f'index 0:{idx} and {tmres}'
        else:
            ch += f'index {idx+1}:{idx*2} and {tmres}'
        ch += ') or ('
    sel = ch[:-4]
    print(sel)
    ## Get the centroid of the whole trajectory, used as representative structure for the pymol renderings
    ## Long to run, the resulting dictionary is: times = {'wt':786.4, 'M654V':412.1}
    t = get_centroid(universe, sel, path, aligned = True, outname=f'centroid_{var}.pdb')
    times[var] = t
'''
vars = ['wt', 'M654V'] # 
times = {'wt':786.4, 'M654V':412.1}
ions = ['chloride', 'potassium']
scripts = f'{basepath}/scripts/'
for var in vars:
    path = f'{basepath}trj_data/{var}/xtcs/'
    if f'centroid_{var}.pdb' not in listdir(path):
        universe = mda.Universe(f'{path}TMC1_{var}_p_l.pdb', f'{path}TMC1_{var}_p_l.xtc')
        myframe = [i.frame for i in universe.trajectory if i.time/1000 == times[var]][0]
        print(myframe)
        universe.trajectory[myframe]
        print(universe.trajectory.time/1000, universe.trajectory.frame)
        ats = universe.select_atoms('protein or resname POPC')
        ats.write(f'{path}centroid_{var}.pdb', bonds="conect")
    ## Make the directory for the output
    hole_csvs_folder = f'{basepath}/hole/csv/{var}/'
    dipoles_csvs_folder = f'{basepath}/dipoles/'
    msd_folder = f'{basepath}msd/{var}/'
    outdir = f'{renderings}/'

    try:
        mkdir(outdir)
    except:
        True

    for ion in ions:
        lowD_file = f'{msd_folder}/{ion}/resulting_dfs/positions_w_lowD.csv'
        threshold = 0.03    
        # This threshold corresponds to approximately:
        #                       7/1000 slowest diffusing K ions in wt
        #                       9/1000 slowest diffusing Cl ions in wt
        #                       1% slowest diffusing K ions in M654V
        #                       1.05% slowest diffusing Cl ions in M654V
        # With a threshold around 0.3 you get instead approx the 5% slowest ions. With 0.03 it is  possible to appreciate the presence of potassium ions 
        # interacting with POPC heads.
        #filtered_df = filter_positions(lowD_file, threshold)
        #filtered_df.to_csv(f'{outdir}/filtered_lowD_{ion}_{var}.csv', header=True, index=False)


    # Get selection of pore-lining atoms and the structure with the VdW centres
    centroid = f'{path}centroid_{var}.pdb'
    get_pymol_representations(var, times[var], hole_csvs_folder, dipoles_csvs_folder, centroid, outdir, scripts, monomers = ['A', 'B'], around_membrane = (70, 130), dist = 100)
