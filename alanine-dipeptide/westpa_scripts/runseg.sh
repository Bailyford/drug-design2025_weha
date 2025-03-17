
#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF
tempF=$(tail -n 1 $WEST_SIM_ROOT/common_files/lambda.dat | awk '{print $2}')
python $WEST_SIM_ROOT/common_files/modprmtop.py $tempF CHARGE $WEST_SIM_ROOT/common_files/diala.prmtop ./mod.prmtop 22
python $WEST_SIM_ROOT/common_files/modprmtop.py $tempF DIHEDRAL_FORCE_CONSTANT ./mod.prmtop ./mod.prmtop
python $WEST_SIM_ROOT/common_files/modprmtop.py $tempF LENNARD_JONES_ACOEF ./mod.prmtop ./mod.prmtop 7 10
python $WEST_SIM_ROOT/common_files/modprmtop.py $tempF LENNARD_JONES_BCOEF ./mod.prmtop ./mod.prmtop 7 10
 
#python $WEST_SIM_ROOT/common_files/modify_data_with_lambda.py $tempF $WEST_SIM_ROOT/common_files/p1_nowat.prmtop ./p1_nowat_mod.prmtop


ln -sv $WEST_SIM_ROOT/common_files/diala.prmtop .
ln -sv $WEST_SIM_ROOT/common_files/mod.prmtop .
ln -sv $WEST_SIM_ROOT/common_files/diala.pdb .

#echo $WEST_PARENT_DATA_REF

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/production.in > production.in
#  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.ncrst ./parent.ncrst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/production.in > production.in
  #sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  cp $WEST_PARENT_DATA_REF.ncrst ./parent.ncrst
  #ln -sv $WEST_PARENT_DATA_REF/bstate.ncrst ./parent.ncrst
fi

export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES


#tempP=$(awk -v "iter=$WEST_CURRENT_ITER" 'NR==iter' $WEST_SIM_ROOT/common_files/temp.dat | awk '{print $2}')


echo $tempF $WEST_SIM_ROOT

$PMEMD  -O -i production.in   -p mod.prmtop -c parent.ncrst \
           -r seg.ncrst -x seg.nc      -o seg.log    -inf seg.nfo -AllowSmallBox

#DIST=$(mktemp)
COMMAND="         parm mod.prmtop\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND autoimage \n"
COMMAND="$COMMAND strip :WAT \n"
COMMAND="$COMMAND energy  @1-22  out energy.dat \n"
COMMAND="$COMMAND trajout nowater.nc \n"
COMMAND="$COMMAND go\n"
python $WEST_SIM_ROOT/common_files/get_dihedrals.py
echo -e $COMMAND | $CPPTRAJ
#python $WEST_SIM_ROOT/common_files/shift_energy.py
cat phi.dat > $WEST_PCOORD_RETURN
cat psi.dat > $WEST_PSI_RETURN
cat energy.dat | tail -n +2 | awk '{print $4}' > $WEST_DIH_ENERGY_RETURN
cat energy.dat | tail -n +2 | awk '{print ($5 + $7)}' > $WEST_VDW_ENERGY_RETURN
cat energy.dat | tail -n +2 | awk '{print ($6 + $8)}' > $WEST_ELEC_ENERGY_RETURN
cat energy.dat | tail -n +2 | awk '{print $2}' > $WEST_BOND_ENERGY_RETURN
cat energy.dat | tail -n +2 | awk '{print $3}' > $WEST_ANGLE_ENERGY_RETURN
python $WEST_SIM_ROOT/common_files/get_energy.py > $WEST_ENERGY_RETURN
#cat $RMSD | tail -n +2 | awk {'print $2 , $3'} > $WEST_AUX_RETURN
python $WEST_SIM_ROOT/common_files/get_coordinates.py
