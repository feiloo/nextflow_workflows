# Modules

ukb_main_workflow implements a main workflow that delegates to one of the subworkflows via the --workflow_variation parameter.

arriba_nextflow allows running the arriba tool in a simple reproducable fashion.

variantinterpretation implements a workflow for analysing variants from vcf files.

sequence_alignment implements an variant_callig workflow from fastq to vcf.

sclust is not fully implemented yet.

## usage

### run arriba

```
nextflow run ukb_main_workflow/  \
	-c ~/nextflow_conf_general.config \
	-with-report /path/nextflow_report.html \
	-with-timeline /path/timeline.html \
	-with-trace /path/trace.txt \
	--reference_data /path/arriba_reference/GRCh38+GENCODE38/ \
	--samplesheet /path/input_Arriba/samplesheets/samples_DD_MM_YYYY.csv \
	--output_dir /path/output_Arriba/workflow_run_DD_MM_YYYY/ 
```

### run sequence_alignment

```
cd $NEXTFLOW_CALLDIR && nextflow run $NEXTFLOW_MODULES/ukb_main_workflow \
	-c $NEXTFLOW_CONFIGS_CUSTOM/nextflow_conf_specific.config \
	--workflow_variation sequence_alignment \
	--samplesheet $PRIVATE_TESTDATA_DIR/samplesheets/samplesheet_wes_ukb_main_workflow.csv \
	-resume
