#!/bin/bash
vars='wt M654V'
basepath='/home/davide/storage/lab_dati/davide_TMC1/definitivi/'
outfolder=$basepath/'rdf'/
scripts=$basepath/'scripts'/
mkdir $outfolder

for var in $vars; do

out=$outfolder/$var/
mkdir $out

traj_folder=$basepath'trj_data/'$var'/xtcs/'
struc=$basepath'trj_data/'$var'/TMC1_'$var'_npt.tpr'
traj=$traj_folder'TMC1_'$var'_noPBC.xtc'
index=$basepath'trj_data/'$var'/index.ndx'
### To calculate Radial Distribution Functions:
baseout=$out'rdf_TMC1_'$var
## Chloride

gmx rdf -f $traj -s $struc -ref Protein -selrpos whole_mol_com -sel POPC -seltype whole_mol_com -n $index -xy -o $baseout'_POPCs_xy'.xvg
done

python3 $scripts/rdf_plots.py $outfolder
