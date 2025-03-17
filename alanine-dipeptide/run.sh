#!/bin/bash
#SBATCH --job-name=diala
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=48:00:00
#SBATCH --cluster=teach
##SBATCH --mail-type=FAIL,COMPLETE
##SBATCH --partition=lchong
##SBATCH --partition=

# run.sh
#
# Run the weighted ensemble simulation. Make sure you ran init.sh first for bash.
#
#set -x
#cd $SLURM_SUBMIT_DIR
#module load gcc/8.2.0 openmpi/4.0.3
#module load amber/22

#source env.sh || exit 1
#env | sort
#cd $WEST_SIM_ROOT
#source init.sh
# BASH, make sure init.sh is executed first
#w_run --work-manager zmq "$@" &> west-${SLURM_JOBID}.log
echo "START `date +"%Y"-%m-%dT%T` `date +"%s"`"
echo "JobID:    ${SLURM_JOB_ID}"
echo "User:     ${USER}"
echo "Hostlist: ${SLURM_NODELIST}"
echo "$CUDA_VISIBLE_DEVICES"
#-------------------------
# WE INIT & RUN
  source env.sh || exit 1
  cd $WEST_SIM_ROOT
  SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json
  source init.sh && echo 'yay'
# start server
  #_run --debug --work-manager=zmq \
  w_run --verbose --work-manager=zmq \
      --n-workers=0 \
      --zmq-mode=master \
      --zmq-write-host-info=$SERVER_INFO \
      --zmq-comm-mode=ipc \
      --zmq-master-heartbeat 3 \
      --zmq-worker-heartbeat 120 \
      --zmq-startup-timeout 3600 \
      --zmq-shutdown-timeout 360 \
      --zmq-timeout-factor 240 \
      &> west-$SLURM_JOBID.log &

  echo " -> Finished launching master."

# wait on master to create SERVER_INFO file
  for ((n=0; n<180; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
  done

  echo " -> Found SERVER_INFO file, launching clients."

# exit if host info file doesn't appear in one minute
  if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
  fi

  echo " -> Found SERVER_INFO file, launching clients."

  scontrol show hostname $SLURM_NODELIST >& SLURM_NODELIST.log

# start clients on nodes, with the proper number of cores on each
  for node in $(scontrol show hostname $SLURM_NODELIST); do
    echo " -> Launching client on $node."
    export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh \
        $SLURM_SUBMIT_DIR \
        $SLURM_JOBID \
        $node \
        $CUDA_VISIBLE_DEVICES \
         --verbose \
          --work-manager=zmq \
          --zmq-mode=client \
          --zmq-read-host-info=$SERVER_INFO \
          --zmq-comm-mode=tcp \
          --n-workers=$SLURM_NTASKS_PER_NODE &> node.log$$ &
  done
  echo " -> Launching Completed... ."

  wait

#-------------------------

echo "END `date +"%Y"-%m-%dT%T` `date +"%s"`"

