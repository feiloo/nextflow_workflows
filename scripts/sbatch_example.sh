#!/bin/bash

#SBATCH --job-name='oncoscanner'
#SBATCH --mem 16000 
#SBATCH --time 240:00:00
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --error="slurm-%j.err"
#SBATCH --output="slurm-%j.out"
#SBATCH --nice=0

set -euo pipefail

export RUNDIR_CUSTOM=$1 
export NXF_CACHE_DIR=$RUNDIR_CUSTOM/nxf_cachedir 
export NEXTFLOW_SCRATCH_DIR=$RUNDIR_CUSTOM/scratch

# usage hint: sbatch run_routine.sh $(pwd)/rundir8 samplesheet.csv md5sum_renamed.txt 'wgs'
# use this instead of nextflow versions (nextflow run -r), because nextflow versioning doesnt work with local git repos
# and we have to use local git repos for availability and because nextflow doesnt support running repos with non toplevel nextflow.configs
COMMIT=8c44c36ff39cf4da79257e693b94277a15815ff9

# ensure the nextflow_workflows git repo is clean, and everything is commited
pushd $NEXTFLOW_MODULES
git diff --quiet && git diff --cached --quiet
git checkout --recurse-submodules $COMMIT
git diff --quiet && git diff --cached --quiet
popd

# copy to local path to ensure the modules arent modified
TARGET_DIR=$RUNDIR_CUSTOM/tmp/

if [[ -d "$TARGET_DIR" ]]; then
	echo "Target directory '$TARGET_DIR' already exists"
	exit 1
fi

mkdir -p $TARGET_DIR
cp --no-clobber -r $NEXTFLOW_MODULES/../ $TARGET_DIR
export NEXTFLOW_MODULES=$TARGET_DIR/modules

SAMPLESHEET=$RUNDIR_CUSTOM/$2
MD5SUMS=$RUNDIR_CUSTOM/$3
LIBRARY_TYPE=$4

export TMPDIR=$RUNDIR_CUSTOM/tmpdir/
export NEXTFLOW_WORKDIR_CUSTOM=$RUNDIR_CUSTOM/nextflow_workdir

echo params are $1 $2 $3 $4
echo path is $(pwd)

# print current sbatch script to out
echo '---- BEGIN SCRIPT ---- '
cat "$0"
echo '---- END SCRIPT ---- '
echo '---- BEGIN env ---- '
printenv
echo '---- END env ---- '
echo '---- BEGIN CONFIG ---- '
cat $NEXTFLOW_MODULES/oncoscanner/slurm.config
echo '---- END CONFIG ---- '

# use exec to propagate kill-signals, see: https://dhruveshp.com/blog/2021/signal-propagation-on-slurm/
exec nextflow \
	-log nextflow-$SLURM_JOB_ID.log \
	run $NEXTFLOW_MODULES/oncoscanner/ \
        --workflow_variation align_interpret \
	--samplesheet $SAMPLESHEET \
	--hash_db $MD5SUMS \
	-c $NEXTFLOW_MODULES/oncoscanner/slurm.config \
        --library_type $LIBRARY_TYPE \
	--input_dir $RUNDIR_CUSTOM \
	-profile standard \
	-offline \
	-resume \
	-dump-hashes json \
	--tag routine_establish,$LIBRARY_TYPE,$COMMIT \
        --bwa_tool bwa2

PIPELINE_EXIT=$?
echo "finished pipeline with code $PIPELINE_EXIT"
echo cache $NXF_CACHE_DIR workdir $NEXTFLOW_WORKDIR_CUSTOM 
