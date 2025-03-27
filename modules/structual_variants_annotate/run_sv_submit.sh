#!/bin/bash
path_sv=input_dir

nextflow run $path_sv/sv_ready_for_submit.nf --input_csv $path_sv/input.csv \
                                             --genelist "/projects/wgs_pilot/sv/Genliste_somatisch.csv" \
                                             --outdir $path_sv/process_out_sv -with-conda -resume
