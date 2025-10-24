#!/bin/bash

#SBATCH --job-name='nextflow_workflow'
#SBATCH --mem 16000 
#SBATCH --time 240:00:00
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --error="slurm-%j.err"
#SBATCH --output="slurm-%j.out"
#SBATCH --nice=0

# export RUNDIR_CUSTOM=$(pwd)
export RUNDIR_CUSTOM=/home/pipeline_user/nextflow_calldir/

# ensure the nextflow_workflows git repo is clean, and everything is commited
pushd $NEXTFLOW_MODULES
git diff --quiet && git diff --cached --quiet
popd

# use exec to propagate kill-signals, see: https://dhruveshp.com/blog/2021/signal-propagation-on-slurm/
exec nextflow-25.04.3 run $NEXTFLOW_MODULES/oncoscanner/ \
        --workflow_variation align_interpret \
	--samplesheet $RUNDIR_CUSTOM/samplesheet.csv \
	--hash_db $RUNDIR_CUSTOM/md5sums.txt \
	-c $NEXTFLOW_MODULES/oncoscanner/user.config \
	-profile standard \
	-offline \
	-resume \
	--tag test,wes \
        --bwa_tool bwa2
