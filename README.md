# TMC
Files used for production and analysis of MD trajectories of TMC1 wt and M654V

Concerning the topology, note that the organization of include topology (.itp) files found in ./trj_data/wt/toppar/ and ./trj_data/M654V/toppar/ is different than the ones found in common forcefileds as a result of the automated routine of CHARMM-GUI v3.7. 
The main differences in the forcefield structure are:
 - the files ffbonded.itp and ffnonbonded.itp are absent as the parameters are directly specified in forcefield.itp
 - there are no position restraint (posre_*.itp) files as these are specified under an if statement in PROA.itp (protein chain A), PROB.itp (protein chain B) and POPC.itp

 - ./trj_data folder:
    - TMC1_original.pdb: starting structure
    - wt / M654V folders: input files obtained from CHARMM-GUI, index files, GROMACS .mdp files containing parameters for minimization, NVT and NPT equilibration and production runs and forcefield parameters in the toppar directory (see above)

 - msd folder: data for evaluation of mean square displacement of ions
 - scripts folder: folder with scripts employed for the subsequent analysis:
    - cgo_arrow.py / cgo_grid.py: scripts to display vectors and grids in PyMol
    - default_name_charges.json: dictionary linking atom name to its partial charge
    - density_plotter.py: employed to plot density maps
    - get_indexes.py / get_densmap_indexes.py: scripts to obtain index file for density map and trajectory aligners
    - xpm2matrix.py: to obtain matrices in simpe format from .xpm GROMACS output
    - rdf_plots.py: to plot radial distribution function
    - get_representative_TMC1_traj.py: to obtain reduced trajectory with 1 frame per ns
    - cleaned_hole_analyis.py: perform HOLE analysis and plot the results
    - densmaps.sh: obtain density maps and plot the results
    - dipoles.py: get and plot dipole data
    - get_pymol_data.py: to obtain .pml scripts for visualiation
    - get_surroundings.py: analysis of the surroundings of residue 654
    -  gmx_msd.sh / msd_vis.py: obtain mean square displacement of ions and plot it
    -  rdf.sh: get radial distribution function
