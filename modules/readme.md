search for fastq in pat sequencer:

find sequencer_output_path --maxdepth 2 -name SampleSheet.csv
find sequencer_output_path  -maxdepth 1 -iname SampleSheet.csv -exec cat {} + | cat > index

create the csv:
/path/input_Arriba/samplesheets/samples_DD_MM_YYYY.csv 
with the header: sampe_id,read1,read2

mkdir /path/output_Arriba/workflow_run_DD_MM_YYYY/

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


