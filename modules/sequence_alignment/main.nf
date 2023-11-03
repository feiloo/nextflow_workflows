nextflow.enable.dsl=2

include { fastqc } from "$NEXTFLOW_MODULES/sequence_alignment/fastqc.nf"
include { fastp } from "$NEXTFLOW_MODULES/sequence_alignment/fastp.nf"
include { bwamem2_index_refgenome } from "$NEXTFLOW_MODULES/sequence_alignment/bwamem2.nf"
include { bwamem2_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwamem2.nf"


process stage_fastq {
    // stage fastq files
    // it renames them to the required naming scheme if necessary
    // it also allows copying, in case the datasource is on a slow network drive

    // stageInMode 'copy'

    input:
    tuple val(sample_id), path(normal_read1), path(normal_read2), path(tumor_read1), path(tumor_read2)

    output:
    tuple val(sample_id), path("${sample_id}_N_1.fq.gz"), path("${sample_id}_N_2.fq.gz"), path("${sample_id}_T_1.fq.gz"), path("${sample_id}_T_2.fq.gz")

    script:
    /*
    """
    if [[ ! -f ${sample_id}_N_1.fq.gz ]]; then
    mv ${normal_read1} ${sample_id}_N_1.fq.gz
    fi

    if [[ ! -f ${sample_id}_N_2.fq.gz ]]; then
    mv ${normal_read2} ${sample_id}_N_2.fq.gz
    fi

    if [[ ! -f ${sample_id}_T_1.fq.gz ]]; then
    mv ${tumor_read1} ${sample_id}_T_1.fq.gz
    fi

    if [[ ! -f ${sample_id}_T_2.fq.gz ]]; then
    mv ${tumor_read2} ${sample_id}_T_2.fq.gz
    fi
    """
    */

    """
    if [[ ! -f ${sample_id}_N_1.fq.gz ]]; then
    ln -s ${normal_read1} ${sample_id}_N_1.fq.gz
    fi

    if [[ ! -f ${sample_id}_N_2.fq.gz ]]; then
    ln -s ${normal_read2} ${sample_id}_N_2.fq.gz
    fi

    if [[ ! -f ${sample_id}_T_1.fq.gz ]]; then
    ln -s ${tumor_read1} ${sample_id}_T_1.fq.gz
    fi

    if [[ ! -f ${sample_id}_T_2.fq.gz ]]; then
    ln -s ${tumor_read2} ${sample_id}_T_2.fq.gz
    fi
    """
}


workflow sequence_alignment {
  take:
    args
  main:
    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'normal_read1', 'normal_read2', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)

    // cast the types of the channel
    sample_pairs = csv_channel.map{
    	row -> 
	[sample_id:row.sample_id, 
	normal_read1:file(row.normal_read1),
	normal_read2:file(row.normal_read2),
	tumor_read1:file(row.tumor_read1),
	tumor_read2:file(row.tumor_read2)
	]
    }

    staged_sample_pairs = stage_fastq(sample_pairs)

    // still check the rows for the naming scheme
    def check_row = { row ->
        def sid = row[0]
	def expected_row = [
		"${sid}",
		"${sid}_N_1.fq.gz",
		"${sid}_N_2.fq.gz",
		"${sid}_T_1.fq.gz",
		"${sid}_T_2.fq.gz"
		]
	def actual_row = [
		"${sid}",
		row[1].getName(),
		row[2].getName(),
		row[3].getName(),
		row[4].getName(),
		]
        return actual_row == expected_row
    }

    
    staged_sample_pairs.subscribe{ row ->
  	if(!check_row(row)) {
	  throw new Exception("row doesnt match naming scheme ${row}")
	}
    }
    
    // a flat channel of [sample_id, read] tuples
    all_reads = staged_sample_pairs.flatMap{
	row -> [[row[0], row[1]],
		[row[0], row[2]],
		[row[0], row[3]],
		[row[0], row[4]],
		]
	}

    //qualities = fastqc(all_reads)

    // a flat channel of [sample_id, read1, read2, file_prefix] tuples
    sample_reads_w_prefix = staged_sample_pairs.flatMap{
	row -> [[row[0], row[1], row[2], "${row[0]}_N"],
		[row[0], row[3], row[4], "${row[0]}_T"]
		]
	}

    // output_file_prefix has to be calculated here, 
    // because it has to be known before starting the process as it it an output thereof
    preprocessed_reads = fastp(sample_reads_w_prefix).preprocessed_reads

    bwa_idx = bwamem2_index_refgenome(args.refgenome)

    // bwa_idx.f0123 was has an f prefixed so, the file ending starts with a non-numeral
    // as a convention
    mapped_reads = bwamem2_align(
    	preprocessed_reads, args.refgenome, 
    	bwa_idx.f0123, bwa_idx.amb, 
	bwa_idx.bwt_2bit_64, bwa_idx.ann, 
	bwa_idx.pac)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}


/*
    // filenames have to end with N_1.fastq.gz
    // where N or T is normal or tumor and the number is the read
    // def pattern = ~".*[NT]_[12]\.(fastq|fq)(\.gz)?"

    // the stricter pattern
    // def pattern = ~".*[NT]_[12]\.fastq\.gz"

    // check that filename matches the naming scheme
    for (filename in all_reads.map{row->row[1]).collect()) {
      if (!(filename =~ pattern).find()){
	throw new Exception("filname doesnt match naming scheme ${filename}")
	}
    }
*/
