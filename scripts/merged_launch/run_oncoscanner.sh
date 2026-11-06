#!/bin/bash

#SBATCH --job-name='oncoscanner'
#SBATCH --mem 16000 
#SBATCH --time 280:00:00
#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --error="slurm-%j.err"
#SBATCH --output="slurm-%j.out"
#SBATCH --nice=0
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_hpc_no_oversub

# #SBATCH --nodelist=ulgs004

set -euo pipefail

export RUNDIR_CUSTOM=$1 
#export NXF_CACHE_DIR=$RUNDIR_CUSTOM/nxf_cachedir 
NXF_CACHE_DIR=$NEXTFLOW_CALLDIR/nxf_cachedir_shared
export NXF_CACHE_DIR
export NEXTFLOW_SCRATCH_DIR=$RUNDIR_CUSTOM/scratch
mkdir -p $NEXTFLOW_SCRATCH_DIR

# futher clc subworkflow vars, stubbed out here:
export CLC_WORKFLOW_NAME='none'
export CLC_IMPORT_DIR='none'
export CLC_EXPORT_DIR='none'
export CLC_DESTDIR='none'

# usage hint: sbatch run_routine.sh $(pwd)/rundir8 samplesheet.csv md5sum_renamed.txt 'wgs'
# use this instead of nextflow versions (nextflow run -r), because nextflow versioning doesnt work with local git repos
# and we have to use local git repos for availability and because nextflow doesnt support running repos with non toplevel nextflow.configs

#COMMIT=0f351ead233c46a296ec63f871ec246d2b3b04ab
#COMMIT=ae2c2c6255cf996046b730bcb613fba86b8fbb43
#COMMIT=d00ba910baaf7d75346871d1f1d90d041df53c23
#COMMIT=bf213c39875f4c60e6ddd62d3e62603f1f6d48d7

#prod
#COMMIT=4146625992768597c44a275db8a8874c4fd515d9

# hyperspeed test
#COMMIT=b3b143292f1e033d3127edb769bf7274397815ea
# add zyto and cfdna tumor modalities
#COMMIT=560a6771846dbcbd3d2026914240cf0845094a16
# rename output folders
# COMMIT=eca7b7e0e3e142cea6a6c04b4654e18900a4c783

# greatly increase UKB_tool performance
COMMIT=e5cb4783771cbe00861c14736274edf955ef9ec3

nextflow secrets set oncokb_token $(cat /etc/oncokb_token.txt)

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
fi

mkdir -p $TARGET_DIR
cp --no-clobber -r $NEXTFLOW_MODULES/../../nextflow_workflows $TARGET_DIR/
export NEXTFLOW_MODULES=$TARGET_DIR/nextflow_workflows/modules


SAMPLESHEET=$RUNDIR_CUSTOM/$2
MD5SUMS=$RUNDIR_CUSTOM/$3
LIBRARY_TYPE=$4

export TMPDIR=$RUNDIR_CUSTOM/tmpdir/
export NEXTFLOW_WORKDIR_CUSTOM=$RUNDIR_CUSTOM/nextflow_workdir
UUID_NAME=$RUNDIR_CUSTOM
export UUID_NAME

# this enables us parallel runs and auto-resume
# very cool and important for setting the nextflow resume and db entries correctly, reflects identity for caching:
# 8f596ff1-1987-4ca5-a953-6f296849b682 is now our single namespace uuid
NAMESPACE_UUID="8f596ff1-1987-4ca5-a953-6f296849b682"
export NAMESPACE_UUID

# derive deterministic session uuid from namespace and run name

SESSION_ID="$(python3 - <<'PY'
import uuid
import os
namespace = uuid.UUID(os.environ['NAMESPACE_UUID'])
name = os.environ['UUID_NAME']
print(uuid.uuid5(namespace, name))
PY
)" 

if [[ -z "${SESSION_ID+x}" || -z "$SESSION_ID" ]]; then
    echo "Error: SESSSION_ID is unset or empty" >&2
    exit 1  # or return 1 if in a function
fi


export SESSION_ID

# set nextflow UUID env var to name the next run
NXF_UUID=$SESSION_ID
export NXF_UUID

# automatically add the resume $SESSION_ID if the session exists in nxf_cache_dir
RESUME_FLAGS=()
[[ -d $NXF_CACHE_DIR/cache/$SESSION_ID ]] && RESUME_FLAGS+=(-resume "$SESSION_ID")

echo SESSION_ID is: "$SESSION_ID"
echo resume is: "${RESUME_FLAGS[@]}"

echo params are $1 $2 $3 $4
echo path is $(pwd)


#CONFIG=$NEXTFLOW_CALLDIR/user.config
#CONFIG=$NEXTFLOW_MODULES/oncoscanner/slurm.config
CONFIG=/etc/oncoscanner/oncoscanner.config

# print current sbatch script to out
echo '---- BEGIN SCRIPT ---- '
cat "$0"
echo '---- END SCRIPT ---- '
echo '---- BEGIN env ---- '
printenv
echo '---- END env ---- '
echo '---- BEGIN CONFIG ---- '
#cat $NEXTFLOW_MODULES/oncoscanner/slurm.config
#cat $NEXTFLOW_CALLDIR/slurm.config
cat $CONFIG
echo '---- END CONFIG ---- '

# use srun or use exec to propagate kill-signals, see: https://dhruveshp.com/blog/2021/signal-propagation-on-slurm/
srun --nice=0 --nodelist=ulgs004 --mem=14G --partition=cpu_hpc_no_oversub --cpus-per-task=1 --tasks=1 --time=16-00:00:00 --export=ALL \
	nextflow \
	-log nextflow-$SLURM_JOB_ID.log \
	run $NEXTFLOW_MODULES/oncoscanner/ \
        --workflow_variation align_interpret \
	--samplesheet $SAMPLESHEET \
	--hash_db $MD5SUMS \
	-c $CONFIG \
        --library_type $LIBRARY_TYPE \
	--input_dir $RUNDIR_CUSTOM \
	-profile standard \
	--calculate_tmb false \
	-offline \
	"${RESUME_FLAGS[@]}" \
	--use_proprietary true \
	--tag routine_establish,$LIBRARY_TYPE,$COMMIT \
        --bwa_tool bwa2

# -resume \
# --bedfile $NEXTFLOW_MODULES/../assets/QIAseqExomehg380720.bed \

PIPELINE_EXIT=$?
echo "finished pipeline with code $PIPELINE_EXIT"
echo cache $NXF_CACHE_DIR workdir $NEXTFLOW_WORKDIR_CUSTOM 
