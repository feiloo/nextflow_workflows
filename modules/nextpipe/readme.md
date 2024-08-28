In this nxf workflow CLC PanCancer data obtained from CLC Genomis Server/Workbench are processed. 

SETUP
        bin folder with PAN_varianten_v1.0.py and folder functions_pan has to be in work dir

STEP1
	Export vcf file according to this example (xxxx.VCF.vcf) in folder X:\PAT-Sequenzer\PANCANCER_INPUT_NXF\vcf_input
	Export txt/csv file according to this example(xxxx.CSV.csv) in folder X:\PAT-Sequenzer\PANCANCER_INPUT_NXF\csv_input

STEP2
	Run/execute the follofing command

nextflow run PANCANCER.nf \
         --vcfs  '/PAT-Sequenzer/PANCANCER_INPUT_NXF/vcf_input/*' \
         --clc_csvs '/PAT-Sequenzer/PANCANCER_INPUT_NXF/csv_input/*' \
         --dir_cache /data2/basitta/vep_cache/ \
         --fasta /data2/basitta/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         --transcript_lst /data2/basitta/pan_data/finale_aktuelle_pancancer_transcript_lst.xlsx \
         --outdir /PAT-Sequenzer/PANCANCER_OUTPUT_NXF/ \
         -with-conda \ 
         -with-report xxx (right now: optional)

