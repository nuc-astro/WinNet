#!/bin/bash
# Authors: M. Ugliano, D. Martin

# Task name
#SBATCH -J My_WinNet_run

# In case of a delayed start:
# #SBATCH --begin=now+2hours

# Include time limit here (24h)
#SBATCH --time=24:00:00

# The partition has to match the cluster
#SBATCH --partition=Compute

# Mail address
# #SBATCH --mail-type=end
# #SBATCH --mail-user=here could be your mail

# Working directory on shared storage
#SBATCH -D /path/to/run

# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log

# Function to print log messages
_log() {
  local format='+%Y/%m/%d-%H:%M:%S'
  echo [`date $format`] "$@"
}


tracers=$(ls -d */)
ntracers=$(echo "$tracers" | wc -l)

i=$((SLURM_ARRAY_TASK_ID))
while [ "$i" -le "$ntracers" ]; do
    #rundir=`ls -d */ | sed -n "$i"p`
    rundir=`echo "$tracers" | sed -n "$i"p`
    ######_log "$rundir"
    # check if output exists
    if [ ! -e "$rundir"/OUT ];  then
      echo "$SLURM_JOB_ID:$SLURM_ARRAY_TASK_ID" >> "$rundir"/blocked
      _log Task $SLURM_ARRAY_TASK_ID \($SLURM_JOB_ID\) trying to run on "$rundir"
      sleep 0.3s
    fi
    a=$(head -n1 "$rundir"/blocked)
    # continue looping if some other run was first
    if [[ "$SLURM_JOB_ID:$SLURM_ARRAY_TASK_ID" != "$a" ]]; then
      i=$((i+1))
      continue
    fi

    cd "$rundir"
    _log Task $SLURM_ARRAY_TASK_ID \($SLURM_JOB_ID\) running in "$rundir" on $SLURM_JOB_NODELIST
    ulimit -s unlimited

    #####################################
    ### define other things to do here: #
    #####################################
    export OMP_NUM_THREADS=1
    # Setting a timeout of 2 hours per run just in case something get stuck
    timeout 2h ./winnet *.par >OUT 2>ERR
    rm winnet
    ###############################

    cd ..

    _log Task $SLURM_ARRAY_TASK_ID \($SLURM_JOB_ID\) finished in "$rundir" on $SLURM_JOB_NODELIST
    i=$((i+1))
done

_log Task $SLURM_ARRAY_TASK_ID \($SLURM_JOB_ID\) done.
