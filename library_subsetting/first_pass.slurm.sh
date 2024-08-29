#!/bin/bash

# run: export DATASET_FILENAME='...'; sbatch /opt/xchem-fragalysis-2/mferla/library_making/first_pass.slurm.sh

#SBATCH --job-name=chunk
#SBATCH --chdir=/opt/xchem-fragalysis-2/mferla
#SBATCH --output=/opt/xchem-fragalysis-2/mferla/library_making/logs/slurm-error_%x_%j.log
#SBATCH --error=/opt/xchem-fragalysis-2/mferla/library_making/logs/slurm-error_%x_%j.log
#SBATCH --partition=main
#SBATCH --export=DATASET_FILENAME
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --priority=-666

# -------------------------------------------------------

export SUBMITTER_HOST=$HOST
export HOST=$( hostname )
export USER=${USER:-$(users)}
export DATA=/opt/xchem-fragalysis-2
export HOME=$DATA/mferla
source /etc/os-release;

echo "Running $SLURM_JOB_NAME ($SLURM_JOB_ID) as $USER in $HOST which runs $PRETTY_NAME submitted from $SUBMITTER_HOST"
echo "Request had cpus=$SLURM_JOB_CPUS_PER_NODE mem=$SLURM_MEM_PER_NODE tasks=$SLURM_NTASKS jobID=$SLURM_JOB_ID partition=$SLURM_JOB_PARTITION jobName=$SLURM_JOB_NAME"
echo "Started at $SLURM_JOB_START_TIME"
echo "job_pid=$SLURM_TASK_PID job_gid=$SLURM_JOB_GID topology_addr=$SLURM_TOPOLOGY_ADDR home=$HOME cwd=$PWD"

export CONDA_PREFIX=$DATA/mferla/waconda-slurm
export CONDA_QUIET=true
export CONDA_YES=true
export SLACK_WEBHOOK='https://hooks.slack.com/services/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'

source $CONDA_PREFIX/etc/profile.d/conda.sh
export CONDA_ENVS_PATH="$CONDA_PREFIX/envs:$DATA/mferla/waconda/envs"

conda activate base;

echo "************************"
echo "HELLO $SLURM_JOB_NAME!"
echo "************************"
echo "Greetings from $SLURM_JOB_NAME script ${0} as $USER in $HOST which runs on $PRETTY_NAME with $CONDA_PREFIX"

# END OF FLUFF
# -------------------------------------------------------
# ACTUAL COMMAND

python /opt/xchem-fragalysis-2/mferla/library_making/first_pass.py $DATASET_FILENAME

# END OF ACTUAL COMMAND
# -------------------------------------------------------
# FINAL FLUFF

curl -X POST -H 'Content-type: application/json' --data '{"text":"'$SLURM_JOB_NAME' complete"}' $SLACK_WEBHOOK