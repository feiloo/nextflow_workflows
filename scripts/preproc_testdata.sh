PROCS="$(nproc --all)"

pushd $TESTDATA_DIR/tmpdir

fasterq-dump -x --size-check off --temp $TESTDATA_DIR/tmpdir --disk-limit 400000000000 --disk-limit-tmp 400000000000 --split-files $TESTDATA_DIR/SRR7890824.sra
fasterq-dump -x --size-check off --temp $TESTDATA_DIR/tmpdir --disk-limit 400000000000 --disk-limit-tmp 400000000000 --split-files $TESTDATA_DIR/SRR7890827.sra

# optionally subsample the data
#parallel "seqtk sample -s100 {} 10000 | bgzip -@ $PROCS > subsampled_{}" ::: $TESTDATA_DIR/*.fastq


bgzip -@ $PROCS --keep $TESTDATA_DIR/tmpdir/SRR7890824_1.fastq
bgzip -@ $PROCS --keep $TESTDATA_DIR/tmpdir/SRR7890824_2.fastq
bgzip -@ $PROCS --keep $TESTDATA_DIR/tmpdir/SRR7890827_1.fastq
bgzip -@ $PROCS --keep $TESTDATA_DIR/tmpdir/SRR7890827_2.fastq

#gzip --keep $TESTDATA_DIR/tmpdir/SRR7890824_1.fastq
#gzip --keep $TESTDATA_DIR/tmpdir/SRR7890824_2.fastq
#gzip --keep $TESTDATA_DIR/tmpdir/SRR7890827_1.fastq
#gzip --keep $TESTDATA_DIR/tmpdir/SRR7890827_2.fastq

# change into naming scheme
mv $TESTDATA_DIR/tmpdir/SRR7890824_1.fastq $TESTDATA_DIR/tmpdir/SRR7890824-24_T_1.fastq
mv $TESTDATA_DIR/tmpdir/SRR7890824_2.fastq $TESTDATA_DIR/tmpdir/SRR7890824-24_T_2.fastq
mv $TESTDATA_DIR/tmpdir/SRR7890827_1.fastq $TESTDATA_DIR/tmpdir/SRR7890824-24_N_1.fastq
mv $TESTDATA_DIR/tmpdir/SRR7890827_2.fastq $TESTDATA_DIR/tmpdir/SRR7890824-24_N_2.fastq

# generate samplesheet
echo "sample_id,normal_read1,normal_read2,tumor_read1,tumor_read2
SRR7890824-24,$TESTDATA_DIR/tmpdir/SRR7890824-24_N_1.fastq,$TESTDATA_DIR/tmpdir/SRR7890824-24_N_2.fastq,$TESTDATA_DIR/tmpdir/SRR7890824-24_T_1.fastq,$TESTDATA_DIR/tmpdir/SRR7890824-24_T_2.fastq" > $TESTDATA_DIR/test_sequence_alignment_samplesheet.csv
