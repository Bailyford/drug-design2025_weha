#!/bin/bash

# Set up simulation environment
source env.sh

# Clean up from previous/ failed runs
rm -rf traj_segs seg_logs istates west.h5 
mkdir   seg_logs traj_segs istates

#bash common_files/setup_bstates.sh 
# Set pointer to bstate and tstate
BSTATE_ARGS="--bstate-file $WEST_SIM_ROOT/bstates/bstates.txt"
#TSTATE_ARGS="--tstate-file $WEST_SIM_ROOT/tstate.file"

echo ${BSTATE_ARGS}
# Run w_init
w_init \
  $BSTATE_ARGS \
  --segs-per-state 1 \
  --work-manager=threads "$@"
