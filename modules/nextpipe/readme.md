In this nxf workflow CLC PanCancer data obtained from CLC Genomis Server/Workbench are processed. 

This organizes pancancer data into nicely readable tables.

example usage:

```
nextflow run PANCANCER.nf \
         --vcfs  '/PANCANCER_INPUT_NXF/vcf_input/*' \
         --clc_csvs '/PANCANCER_INPUT_NXF/csv_input/*' \
         --dir_cache /data/reference/vep_cache/ \
         --fasta /data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         --transcript_lst transcriptlist.xlsx \
         --outdir /PANCANCER_OUTPUT_NXF/ \
         -with-conda \ 
         -with-report report.html

```
