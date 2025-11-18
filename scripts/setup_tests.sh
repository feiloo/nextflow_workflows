#!/bin/bash
set -euo pipefail

TESTDIR=$(pwd)/tests/
mkdir $TESTDIR
cd $TESTDIR

git clone --recursive git@github.com:nh13/DWGSIM.git

cd DWGSIM
sed -i "s/-lcurses/-lncurses/g" samtools/Makefile

CC=gcc-14 CXX=g++-14 make -j
CC=gcc-14 CXX=g++-14 make test

REFGENOME=$NGS_REFERENCE_DIR/oncoscanner_reference/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

./dwgsim -C 0.01 -r 0.001 $REFGENOME ../simulated_reads_tumor &
P1=$!
./dwgsim -C 0.01 -r 0.0 $REFGENOME ../simulated_reads_normal &
P2=$!

wait $P1
wait $P2

cd ..
mv -v --no-clobber simulated_reads_normal.bwa.read1.fastq.gz X123-25_N_FFPE_1.fq.gz
mv -v --no-clobber simulated_reads_normal.bwa.read2.fastq.gz X123-25_N_FFPE_2.fq.gz
mv -v --no-clobber simulated_reads_tumor.bwa.read1.fastq.gz X123-25_T_FFPE_1.fq.gz
mv -v --no-clobber simulated_reads_tumor.bwa.read2.fastq.gz X123-25_T_FFPE_2.fq.gz

parallel md5sum {} >> md5sum.txt ::: *.fq.gz

echo 'sample_id,normal_modality,tumor_modality' >> samplesheet.csv
echo 'X123-25,FFPE,FFPE' >> samplesheet.csv

