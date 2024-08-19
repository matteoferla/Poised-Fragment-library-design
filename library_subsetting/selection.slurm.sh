#!/bin/bash

# The submitter node is cepheus-slurm.diamond.ac.uk (cs05r-sc-cloud-30.diamond.ac.uk)
# cs05r-sc-cloud-29.diamond.ac.uk is for slurm testing
# wilson is something else?
# submit via `export VC='???'; sbatch /opt/xchem-fragalysis-2/mferla/library_making/sele2.slurm.sh`

#SBATCH --job-name=selection
#SBATCH --chdir=/opt/xchem-fragalysis-2/mferla
#SBATCH --output=/opt/xchem-fragalysis-2/mferla/logs/slurm-error_%x_%j.log
#SBATCH --error=/opt/xchem-fragalysis-2/mferla/logs/slurm-error_%x_%j.log
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --exclusive
##SBATCH --time=24:00:00
#SBATCH --export=DATASET
#SBATCH --nice=100
##SBATCH --export=NONE

# -------------------------------------------------------

export SUBMITTER_HOST=$HOST
export HOST=$( hostname )
export USER=${USER:-$(users)}
export DATA=/opt/xchem-fragalysis-2
export HOME=$DATA/mferla
source /etc/os-release;
export SLACK_WEBHOOK='https://hooks.slack.com/services/👾👾/👾👾/👾👾'


echo "Running $SLURM_JOB_NAME ($SLURM_JOB_ID) as $USER in $HOST which runs $PRETTY_NAME submitted from $SUBMITTER_HOST"
echo "Request had cpus=$SLURM_JOB_CPUS_PER_NODE mem=$SLURM_MEM_PER_NODE tasks=$SLURM_NTASKS jobID=$SLURM_JOB_ID partition=$SLURM_JOB_PARTITION jobName=$SLURM_JOB_NAME"
echo "Started at $SLURM_JOB_START_TIME"
echo "job_pid=$SLURM_TASK_PID job_gid=$SLURM_JOB_GID topology_addr=$SLURM_TOPOLOGY_ADDR home=$HOME cwd=$PWD"

# -------------------------------------------------------

source $DATA/mferla/waconda-slurm/etc/profile.d/conda.sh
#export CONDA_ENVS_PATH="$DATA/mferla/waconda-slurm/envs:$DATA/mferla/waconda/envs"
#conda activate compchem;
conda activate torch;
cd $DATA/mferla/library_making
pwd;
echo $DATASET;

mv $DATASET /tmp/$DATASET
python selection_v2.py /tmp/$DATASET common_synthons_SpikeIn.pkl.gz


curl -X POST -H 'Content-type: application/json' --data '{"text":"selection '$DATASET' complete"}' $SLACK_WEBHOOK



# -------------------------------------------------------

echo 'complete'

exit 0
