#!/bin/bash
set -euo pipefail
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_T_1.bwa.dedup.bam
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_N_1.bwa.dedup.bam


# convert to vastq:
samtools index -@8 WES_EA_T_1.bwa.dedup.bam
samtools index -@8 WES_EA_N_1.bwa.dedup.bam

# # warning, the bams are made from deduped fastqs
# samtools fastq -@8 -n -1 WES_EA_T_1_R1.fastq -2 WES_EA_T_1_R2.fastq \
# 	-0 /dev/null -s /dev/null WES_EA_T_1.bwa.dedup.bam
# 
# # warning, the bams are made from deduped fastqs
# samtools fastq -@8 -n -1 WES_EA_N_1_R1.fastq -2 WES_EA_N_1_R2.fastq \
# 	-0 /dev/null -s /dev/null WES_EA_N_1.bwa.dedup.bam
# 
# bgzip -@ 8 WES_EA_T_1_R1.fastq
# bgzip -@ 8 WES_EA_T_1_R2.fastq
# bgzip -@ 8 WES_EA_N_1_R1.fastq
# bgzip -@ 8 WES_EA_N_1_R2.fastq
