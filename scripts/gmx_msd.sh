bp='/home/davide/storage/lab_dati/davide_TMC1/definitivi/'
outfolder=$bp'msd/'
mkdir $outfolder
protein='TMC1'
vars='wt M654V' #
scripts=$bp'/scripts/'

for var in $vars
do
traj_path=$bp'/trj_data/'$var'/'
ndx=$traj_path'index.ndx' 
trj=$traj_path'TMC1_'$var'_npt.trr'
tpr=$traj_path'TMC1_'$var'_npt.tpr'
out=$outfolder/$var
mkdir $out

### BY TIME: considering just ions
step=100 # ns
tot=1000

o=$out
op=$o/'potassium'/
mkdir $op
oc=$o/'chloride'/
mkdir $oc

let i=$step;
while [ $i -le $tot ]
do
let start1=$i-$step

covgr="14\n"
loc='potassium'
gmx msd -f $trj -n $ndx -s $tpr -o $op/'msdout_'$loc'.xvg' -mol $op/'diff_mol_'$loc'_'$start1'to'$i'.xvg' -b $start1 -e $i -tu ns <<< $covgr

covgr="15\n"
loc='chloride'
gmx msd -f $trj -n $ndx -s $tpr -o $oc/'msdout_'$loc'.xvg' -mol $oc/'diff_mol_'$loc'_'$start1'to'$i'.xvg' -b $start1 -e $i -tu ns <<< $covgr

let i=i+$step
done
done
