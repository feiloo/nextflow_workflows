#!/bin/bash
set -euo pipefail

RUNDIR=$(pwd)/$1
LIBRARY_TYPE=$2

sbatch run_oncoscanner.sh $RUNDIR/oncoscanner samplesheet_oncoscanner.csv md5sum_renamed_short.txt $LIBRARY_TYPE
sbatch run_oa.sh $RUNDIR/oncoanalyser samplesheet_oncoanalyser.csv $LIBRARY_TYPE
