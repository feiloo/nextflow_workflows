#!/bin/bash

nextflow run PANCANCER.nf \
         -c resources.config \
         --vcfs  '/PAT-Sequenzer/PANCANCER/INPUT/vcf_input/*.vcf' \
         --clc_csvs '/PAT-Sequenzer/PANCANCER/INPUT/csv_input/*.csv' \
         --dir_cache /data/reference/vep_cache/111/ \
         --fasta /data/reference/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         --transcript_lst /data2/basitta/pan_data/in_use_finale_aktuelle_pancancer_transcript_lst.xlsx \
         --variantDBi /PAT-Sequenzer/PanCancer_test/Variantenliste22_12_15.xlsx \
         --outdir /PAT-Sequenzer/PANCANCER/OUTPUT/CLC24 \
         -with-conda \
         -with-report /PAT-Sequenzer/PANCANCER/ARCHIVE_AND_LOGS/NXF_REPORTS/$(date +"%Y%m%d%H%M"_nxf.html) \
         -resume

