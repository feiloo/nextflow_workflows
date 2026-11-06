#!/bin/bash
#SBATCH --job-name=oncoanalyser
#SBATCH --mem=16G
#SBATCH --export=ALL
#SBATCH --error="slurm-%j.err"
#SBATCH --output="slurm-%j.out"
#SBATCH --nice=0
#SBATCH --time=17-00:00:00
#SBATCH --partition=cpu_hpc_no_oversub
#SBATCH --nodelist=ulgs004

set -euo pipefail

samplesheet=""
mode=""

show_usage() {
    echo "Usage: $0 -s|--samplesheet <FILE> -m|--mode <MODE>"
    echo "Options:"
    echo "  -s, --samplesheet <FILE>   Path to the samplesheet file (required)"
    echo "  -m, --mode <MODE>          Processing mode: 'wgs', 'wes', (required)"
    echo "  -h, --help                 Show this help message"
    exit 1
}

# Check for no arguments
if [[ $# -eq 0 ]]; then
    show_usage
fi

RUNDIR_CUSTOM=$1
export RUNDIR_CUSTOM
SAMPLESHEET=$RUNDIR_CUSTOM/$2
MODE=$3

NXF_CACHE_DIR=$NEXTFLOW_CALLDIR/nxf_cachedir_shared
export NXF_CACHE_DIR
# export NEXTFLOW_SCRATCH_DIR=$RUNDIR_CUSTOM/scratch
# mkdir -p $NEXTFLOW_SCRATCH_DIR

# stage refdata like: ./nextflow-25.04.7 run nf-core/oncoanalyser   -revision 2.2.0   -profile conda   --mode prepare_reference   --ref_data_types wgs   --genome GRCh38_hmf   --outdir output/

echo rundir $RUNDIR_CUSTOM samplesheet $SAMPLESHEET mode $MODE

# Validate required arguments
if [[ -z "$SAMPLESHEET" ]]; then
    echo "Error: --samplesheet is required" >&2
    exit 1
fi

if [[ -z "$MODE" ]]; then
    echo "Error: --mode is required" >&2
    exit 1
fi

if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "Error: Samplesheet file not found: $samplesheet" >&2
    exit 1
fi

# Validate mode options
case "$MODE" in
    wes|wgs)
        echo "Running in $mode mode using samplesheet: $samplesheet"
        ;;
    *)
        echo "Error: Invalid mode '$mode'. Valid options: wgs, wes" >&2
        exit 1
        ;;
esac

cd $RUNDIR_CUSTOM

# set in case its run as bash script without sbatch
SLURM_JOB_ID="${SLURM_JOB_ID:=-1}"

pseudo_java_instant_time=$(date --utc +"%Y-%m-%dT%H_%M_%S.%NZ")

#pseudo_java_instant_time="2026-04-14T13_53_10.925984586Z"
timestamp="$pseudo_java_instant_time"
RUNNAME="$pseudo_java_instant_time"
echo RUNNAME is $RUNNAME
echo SLURM_JOB_ID is $SLURM_JOB_ID


## cache handling
export TMPDIR=$RUNDIR_CUSTOM/tmpdir/
export NEXTFLOW_WORKDIR_CUSTOM=$RUNDIR_CUSTOM/nextflow_workdir
UUID_NAME="$RUNDIR_CUSTOM"_oncoanalyser
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

echo params are $1 $2 $3
echo path is $(pwd)

# set variables
# just use system default

# set java initial and max heap space
export NXF_JVM_ARGS="-Xms2g -Xmx6g"

#export LANG=en_US.UTF-8
#export LC_ALL=en_US.UTF-8

outdir=$(pwd)/outdirs/$RUNNAME

ONCOANALYSER_REFERENCE_DIR=$NGS_REFERENCE_DIR/oncoanalyser_reference/refdata_copy
export ONCOANALYSER_REFERENCE_DIR=$NGS_REFERENCE_DIR/oncoanalyser_reference/refdata_copy

mkdir -p $NEXTFLOW_WORKDIR_CUSTOM
mkdir -p $outdir

echo workdir is $NEXTFLOW_WORKDIR_CUSTOM
echo outdir is $outdir

echo starting oncoanalyser
echo jobid: $SLURM_JOB_ID

# MODE=wes
# SAMPLESHEET=$samplesheet #$(pwd)/samplesheet_wes.csv

# must be set via outside here bashrc 
# ONCOANALYSER_CONFIG_DIR
#MODE=wgs

if [[ "$MODE" == "wes" ]]; then
	CONFIG=/etc/oncoscanner/oncoanalyser-2.2.0_custom.config
	srun --nice=0 --nodelist=ulgs004 --mem=14G --partition=cpu_hpc_no_oversub --cpus-per-task=1 --tasks=1 --time=16-00:00:00 --export=ALL \
	  nextflow run nf-core/oncoanalyser \
	  -profile conda \
	  -revision 2.2.0 \
	  --mode targeted \
	  --panel wes_mtn \
	  --genome GRCh38_hmf \
	  --input $SAMPLESHEET \
	  --outdir $outdir \
	  --processes_exclude virusinterpreter \
	  -work-dir $NEXTFLOW_WORKDIR_CUSTOM \
	  -dump-hashes \
	  --ref_data_hla_slice_bed $ONCOANALYSER_REFERENCE_DIR/grch38_alt.plus_homologous.bed \
	  -with-conda \
	  -config $CONFIG \
	  "${RESUME_FLAGS[@]}" \
	  --force_panel


# #writes its outputs into the same outputdir additionally


elif [[ $MODE == "wgs" ]]; then

#SBATCH --time=7-00:00:00
	CONFIG=/etc/oncoscanner/oncoanalyser-refdata-2.2.0-wgts.local.config
	srun --nice=0 --nodelist=ulgs004 --mem=14G --partition=cpu_hpc_no_oversub --cpus-per-task=1 --tasks=1 --time=16-00:00:00 --export=ALL \
	  nextflow run nf-core/oncoanalyser \
	  -profile conda \
	  -revision 2.2.0 \
	  --mode wgts \
	  --genome GRCh38_hmf \
	  --input $SAMPLESHEET \
	  --processes_exclude virusinterpreter \
	  --outdir $outdir \
	  -work-dir $NEXTFLOW_WORKDIR_CUSTOM \
	  --ref_data_hla_slice_bed $ONCOANALYSER_REFERENCE_DIR/grch38_alt.plus_homologous.bed \
	  -with-conda \
	  "${RESUME_FLAGS[@]}" \
	  -config $CONFIG

fi

echo successfully finished oncoanalyser, mode: $MODE
echo jobid: $SLURM_JOB_ID

echo starting postprocessing
echo workdir is $NEXTFLOW_WORKDIR_CUSTOM
echo outdir is $outdir

# refdir=oa_refdata
refdir=$NGS_REFERENCE_DIR/oncoanalyser_reference/oa_refdata
# workdir=run4/workdir
# outdir=run4/run_26_09_16

# TODO (only wes mode was tested for now)

cd $NEXTFLOW_CALLDIR
nextflow run postprocessing_wxs_custom_OA.nf \
              -c "/etc/oncoscanner/postprocessing_wxs_custom_OA.config" \
              --wes_bedfile "$refdir/QIAseqExomehg380720.bed" \
              --outdir "$outdir" \
              --ref_genome_fasta "$refdir/ref_oncoanalyser_v2/HMFtools-Resources_ref_genome_38_GRCh38_masked_exclusions_alts_hlas.fasta" \
              --ref_genome_fai "$refdir/ref_oncoanalyser_v2/HMFtools-Resources_ref_genome_38_GRCh38_masked_exclusions_alts_hlas.fasta.fai" \
              --ref_genome_dict "$refdir/ref_oncoanalyser_v2/HMFtools-Resources_ref_genome_38_GRCh38_masked_exclusions_alts_hlas.fasta.dict" \
              --grz_bedfile "$refdir/hg38_440_omim_genes.bed" \
              --wes_allele_files "$refdir/ref_ascat/wes/G1000_alleles_hg38" \
              --wes_loci_files "$refdir/ref_ascat/wes/G1000_loci_hg38" \
              --wes_gc_file "$refdir/ref_ascat/wes/GC_G1000_hg38.txt" \
              --wes_rt_file "$refdir/ref_ascat/wes/RT_G1000_hg38.txt" \
              -work-dir $NEXTFLOW_WORKDIR_CUSTOM \
              --target_dir "$outdir" \
              --library_type $MODE \
              -with-conda

	      # -resume 789cb943-039d-5023-ba5e-63b76b956a32

echo finished postprocessing, jobid: $SLURM_JOB_ID
