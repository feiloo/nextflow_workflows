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
REPO_URL='https://github.com/feiloo/DWGSIM.git'
COMMIT=7ed576b025395c7a376e5f274a6370c529a7500e

pushd "$TMP_DIR" > /dev/null
git clone --branch fast --recurse-submodules "$REPO_URL" "$TMP_DIR/repo"
pushd "$TMP_DIR/repo" > /dev/null
git checkout "$COMMIT"
sed -i 's/-lcurses/-lncurses/g' samtools/Makefile
make -j$(nproc)
popd > /dev/null

ln -s "$(realpath "$REF_GENOME")" "$REF_BASENAME"
samtools faidx "$REF_BASENAME"

# Generate matched tumor‑normal reads (prefix = ${SAMPLE_ID}-25)
"$TMP_DIR/repo/dwgsim" --matched -N 1000000 -z 13 \
    -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
    "$TMP_DIR/$REF_BASENAME" "$DATA_OUTPUTDIR/${SAMPLE_ID}-25"

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
${SAMPLE_ID}-25,BLOOD,FFPE
EOF

# Generate md5sum.txt for the fastq files
echo "Generating md5sum.txt"
cd "$DATA_OUTPUTDIR"
md5sum *.fq.gz > md5sum.txt
cd - > /dev/null

echo "Generated samplesheet.csv, md5sum.txt and fastq files in $DATA_OUTPUTDIR"

