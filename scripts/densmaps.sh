basepath='/home/davide/storage/lab_dati/davide_TMC1/definitivi/'
outpath=$basepath'densmaps/'
mkdir $outpath

protein='TMC1'
vars='wt M654V'
scripts='/home/davide/storage/lab_dati/davide_TMC1/definitivi/scripts/'

## for densmaps of lipids, protein, ions
for var in $vars
do
  traj_path=$basepath'trj_data/'$var'/xtcs/'
  myindex=$traj_path'densmaps_index.ndx'
  structure=$traj_path'TMC1_'$var'_p_l.gro'
  trj=$traj_path'TMC1_'$var'_p_l.xtc'
  if [ ! -f $myindex ]; then
    echo 'Preparing index file'
    /usr/bin/python3 $scripts'get_indexes.py' $structure $trj $traj_path
  else
    echo 'Index file already present!'
  fi
  
  out=$outpath/$var/
  mkdir $out
  k=0
  zones='S4 S5 S6 S7 no_channel POT CLA Lipids protein'
  zone_array=($zones)  # Split the string into an array
  m=${#zone_array[@]}
  while [ $k -lt $m ]
  do
    # Get the k-th element from the array
    loc="${zone_array[$k]}"
    covgr=$k' \n '
    
    axes='x y z'
    for ax in $axes
    do
      outf=$out/$ax'_axis'/
      mkdir $outf
      m_str=$outf$protein'_'$var'_'$loc'_densmap_'$ax.dat
      
      #gmx densmap -f $trj -n $myindex -aver $ax -od $m_str -bin 0.02 <<< $covgr
      rm $outf/\#*
    done
  ((k++))
  done

  
done

img_path=$outpath'/images/'
mkdir $img_path

/usr/bin/python3 $scripts/density_plotter.py $outpath $img_path

## for densmaps of water
for var in $vars
do
  traj_path=$basepath'trj_data/'$var'/xtcs/'
  myindex=$traj_path'densmaps_index_water.ndx'
  structure=$traj_path'TMC1_'$var'_noPBC.gro'
  trj=$traj_path'TMC1_'$var'_noPBC_reduced.xtc'
  if [ ! -f $myindex ]; then
    echo 'Preparing index file'
    /usr/bin/python3 $scripts'get_indexes.py' $structure $trj $traj_path "T"
  else
    echo 'Index file already present!'
  fi
  
  out=$outpath/$var/
  mkdir $out
  k=0
  zones='S4 S5 S6 S7 no_channel POT CLA Lipids protein water'
  zone_array=($zones)  # Split the string into an array
  m=${#zone_array[@]}
  while [ $k -lt $m ]
  do
    # Get the k-th element from the array
    loc="${zone_array[$k]}"
    covgr=$k' \n '
    
    axes='x y z'
    for ax in $axes
    do
      outf=$out/$ax'_axis'/
      mkdir $outf
      m_str=$outf$protein'_'$var'_'$loc'_densmap_'$ax.dat
      
      #gmx densmap -f $trj -n $myindex -aver $ax -od $m_str -bin 0.02 <<< $covgr
      rm $outf/\#*
    done
  ((k++))
  done

  
done

img_path=$outpath'/images/'
mkdir $img_path

/usr/bin/python3 $scripts/density_plotter.py $outpath $img_path
