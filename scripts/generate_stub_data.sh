#!/bin/bash
set -euo pipefail

REF_GENOME="$1"
DATA_OUTPUTDIR="$2"
SAMPLE_ID="$3"
TMP_DIR="${4:-}"

if [ -z "$TMP_DIR" ]; then
    TMP_DIR=$(mktemp -d)
    cleanup() { rm -rf "$TMP_DIR"; }
    trap cleanup EXIT
else
    mkdir -p "$TMP_DIR"
fi

mkdir -p "$DATA_OUTPUTDIR"

REF_BASENAME="$(basename "$REF_GENOME")"
REF_NEWPATH="$DATA_OUTPUTDIR"/"$REF_BASENAME"

ln -s "$(realpath "$REF_GENOME")" $REF_NEWPATH || true
#echo $(pwd)
pushd ./third_party/DWGSIM/samtools
make 
popd
./third_party/DWGSIM/samtools/samtools faidx $REF_NEWPATH


# Generate matched tumor‑normal reads (prefix = ${SAMPLE_ID}-25)
"./third_party/DWGSIM/dwgsim" --matched -N 1000000 -z 13 \
    -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
    $REF_NEWPATH "$DATA_OUTPUTDIR/${SAMPLE_ID}-25"

pushd $DATA_OUTPUTDIR

# Rename output files to desired pattern
for f in "$DATA_OUTPUTDIR"/"${SAMPLE_ID}"-25.*.fastq.gz; do
    base=$(basename "$f")
    newname=$(echo "$base" | sed -E 's/\.(normal|tumor)\.bwa\.read(1|2)\.fastq\.gz/_\1_FFPE_\2.fq.gz/' | sed 's/normal/N/; s/tumor/T/')
    mv "$f" "$DATA_OUTPUTDIR/$newname"
done

# Create samplesheet.csv
cat > "$DATA_OUTPUTDIR/samplesheet.csv" << EOF
sample_id,normal_modality,tumor_modality
${SAMPLE_ID}-25,FFPE,FFPE
EOF

# Generate md5sum.txt for the fastq files
echo "Generating md5sum.txt"
md5sum *.fq.gz > md5sum.txt
popd

echo "Generated samplesheet.csv, md5sum.txt and fastq files in $DATA_OUTPUTDIR"

