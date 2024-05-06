nextflow.enable.dsl=2

include { fastqc } from "$NEXTFLOW_MODULES/sequence_alignment/fastqc.nf"
include { fastp } from "$NEXTFLOW_MODULES/sequence_alignment/fastp.nf"

include { bwamem2_index_refgenome; bwamem2_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwamem2.nf"
include { bwa_index_refgenome; bwa_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwa.nf"

include { bam_stats; bam_depth; sam_to_bam; sort_bam; index_bam; index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary; gatk_markduplicates; gatk_set_tags; gatk_baserecalibrator; gatk_apply_bqsr } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

include { variant_call } from "$NEXTFLOW_MODULES/sequence_alignment/gatk_mutect.nf"

process stage_fastq {
    // stage fastq files
    // it renames them to the required naming scheme if necessary
    // it also allows copying, in case the datasource is on a slow network drive

    //stageInMode 'copy'
    stageInMode 'link'
    cache 'lenient'

    input:
    tuple val(sample_id), path(normal_read1), path(normal_read2), path(tumor_read1), path(tumor_read2)

    output:
    tuple val(sample_id), path(normal_read1), path(normal_read2), path(tumor_read1), path(tumor_read2)

    script:
    assert normal_read1.exists()
    assert normal_read2.exists()
    assert tumor_read1.exists()
    assert tumor_read2.exists()
    """
    echo hello
    #touch ${normal_read1}
    #touch ${normal_read2}
    #touch ${tumor_read1}
    #touch ${tumor_read2}
    """
}

workflow sequence_alignment {
  take:
    args
  main:
    // set default to cleanup process inputs
    // this breaks resuming the workflow, but saves alot of storage
    if(!args.containsKey('cleanup_intermediate_files') || !args['cleanup_intermediate_files']){
    	args['cleanup_intermediate_files'] = false
	}

    if(!args.containsKey('bwa_tool') || !args['bwa_tool']){
    	args['bwa_tool'] = 'bwa2'
	}

    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'normal_read1', 'normal_read2', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: header, skip: 1)

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

    // explicitely stage inputs
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

    // low_mem uses more memory efficient tools
    // for example bwa-mem instead of bwamem2
    // this trades memory for performance here

    if(args.bwa_tool == 'bwa'){
    	bwa_idx = bwa_index_refgenome(args.refgenome)
    } else if(args.bwa_tool == 'bwa2'){
	bwa_idx = bwamem2_index_refgenome(args.refgenome)
    } else {
	throw new Exception("unknown bwa_tool ${args.bwa_tool}")
    }

    refgenome_index = index_fasta(args.refgenome).fasta_index
    refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict

    preprocessed_reads_no_id = preprocessed_reads.map{it -> [it[1], it[2]]} 

    if(args.bwa_tool == 'bwa'){
    	bams = bwa_align(preprocessed_reads_no_id, args.refgenome, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt, bwa_idx.pac, bwa_idx.sa, args.cleanup_intermediate_files)
    } else if(args.bwa_tool == 'bwa2'){
        // bwa_idx.f0123 was has an f prefixed so, the file ending starts with a non-numeral
        bams = bwamem2_align(preprocessed_reads_no_id, args.refgenome, bwa_idx.f0123, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt_2bit_64, bwa_idx.pac, args.cleanup_intermediate_files)
    	//bams = sam_to_bam(sams, args.cleanup_intermediate_files)
    } else {
	throw new Exception("unknown bwa_tool ${args.bwa_tool}")
    }

    marked_bams = gatk_markduplicates(bams).marked_bams
    tagged_bams = gatk_set_tags(marked_bams, args.refgenome, args.cleanup_intermediate_files).tagged_bams

    known_sites_index = gatk_indexfeaturefile(args.known_sites).known_sites_index
    bam_recal_data = gatk_baserecalibrator(tagged_bams, args.refgenome, refgenome_index, refgenome_dict, args.known_sites, known_sites_index)

    // rebuild key for matching with recal-data
    keyed_bams = tagged_bams.map{ it -> ["${it.getSimpleName()}_recal_data", it] }
    keyed_recal_data = bam_recal_data.map{ it -> ["${it.getSimpleName()}", it] }

    bam_w_recal_data = keyed_bams.join(keyed_recal_data).map{ it -> [it[1], it[2]] }

    bam_recalibrated = gatk_apply_bqsr(bam_w_recal_data, args.refgenome, refgenome_index, refgenome_dict)
    bam_w_depth = bam_depth(bam_recalibrated, args.refgenome)
    bams_w_stats = bam_stats(bam_recalibrated, args.refgenome)

    vcfs = variant_call(
	    bam_recalibrated,
	    args.intervals,
	    args.panel_of_normals,
	    args.germline_resource,
	    args.refgenome,
	    args,
	    )

  emit:
    bam = bam_recalibrated
    vcf = vcfs
    bam_depth = bam_w_depth
    bam_stats = bams_w_stats
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
