PROCS="$(nproc --all)"

parallel "seqtk sample -s100 {} 10000 | bgzip -@ $PROCS > subsampled_{}" ::: *.fastq
