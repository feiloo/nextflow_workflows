#!/bin/bash
# Short script to obtain only PASS variants from mutect2 output
FILE_PATH=*.vcf
for file in $FILE_PATH;
do
bcftools view -h $file > filtered_"${file##*/}"
bcftools view -H -f PASS $file >> filtered_"${file##*/}"
done
