# Modules

the module ukb_main_workflow delegates the inputs to one of the subworkflows via the --workflow_variation parameter

sclust only implements the first steps

sarek_wrapper is not functional yet

arriba_nextflow allows running the arriba tool in a simple reproducable fashion

variantinterpretation implements a workflow for analysing variants from vcf files
invoke nextflow 

### how to run arriba
nextflow run ukb_main_workflow/  \
	-c ~/nextflow_conf_general.config \
	-with-report /path/nextflow_report.html \
	-with-timeline /path/timeline.html \
	-with-trace /path/trace.txt \
	--reference_data /path/arriba_reference/GRCh38+GENCODE38/ \
	--samplesheet /path/input_Arriba/samplesheets/samples_DD_MM_YYYY.csv \
	--output_dir /path/output_Arriba/workflow_run_DD_MM_YYYY/ 

### how to run sequence_alignment

cd $NEXTFLOW_CALLDIR && nextflow run $NEXTFLOW_MODULES/ukb_main_workflow \
	-c $NEXTFLOW_CONFIGS_CUSTOM/nextflow_conf_specific.config \
	--workflow_variation sequence_alignment \
	--samplesheet $PRIVATE_TESTDATA_DIR/samplesheets/samplesheet_wes_ukb_main_workflow.csv \
	-resume

### prerequisite system packages

its required that some system packages are installed:

gsl (gnu-scientific library)

running with conda requires a conda installation (mamba and micromamba are not recommended)
you need to add the proper bioconda channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
