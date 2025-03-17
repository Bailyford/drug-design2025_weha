#!/bin/bash
set -x
grep $(basename $WEST_STRUCT_DATA_REF.ncrst) $WEST_SIM_ROOT/bstates/pcoord.init | awk '{print $2}' >  $WEST_PCOORD_RETURN

grep $(basename $WEST_STRUCT_DATA_REF.ncrst) $WEST_SIM_ROOT/bstates/pcoord.init | awk '{print $3}' > $WEST_ENERGY_RETURN

cp $WEST_SIM_ROOT/common_files/diala.prmtop $WEST_TRAJECTORY_RETURN
cp $WEST_STRUCT_DATA_REF.xml $WEST_TRAJECTORY_RETURN

cp $WEST_SIM_ROOT/common_files/diala.prmtop $WEST_RESTART_RETURN
cp $WEST_STRUCT_DATA_REF.ncrst $WEST_RESTART_RETURN/parent.ncrst
echo $(grep $WEST_STRUCT_DATA_REF.ncrst $WEST_SIM_ROOT/bstates/pcoord.init | awk '{print $2}') 

