nextflow.enable.dsl=2

//include { fastqc } from "$NEXTFLOW_MODULES/sequence_alignment/fastqc.nf"
include { fastp } from "$NEXTFLOW_MODULES/sequence_alignment/fastp.nf"

include { bwamem2_index_refgenome; bwamem2_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwamem2.nf"
include { bwa_index_refgenome; bwa_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwa.nf"

include { bam_coverage; bam_stats; bam_depth; sam_to_bam; sort_bam; index_bam; index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary; gatk_markduplicates; gatk_set_tags; gatk_baserecalibrator; gatk_apply_bqsr } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

include { variant_call } from "$NEXTFLOW_MODULES/sequence_alignment/gatk_mutect.nf"

process stage_fastq {
    // stage fastq files
    // it renames them to the required naming scheme if necessary
    // it also allows copying, in case the datasource is on a slow network drive

    //stageInMode 'copy'
    //stageInMode 'link'
    cache 'lenient'

    input:
    tuple val(sample_id), path(normal_read1), path(normal_read2), path(tumor_read1), path(tumor_read2), val(normal_read1_hash), val(normal_read2_hash),  val(tumor_read1_hash),  val(tumor_read2_hash)

    output:
    tuple val(sample_id), path(normal_read1), path(normal_read2), path(tumor_read1), path(tumor_read2), val(normal_read1_hash), val(normal_read2_hash),  val(tumor_read1_hash),  val(tumor_read2_hash)

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

def hash_db_to_dict( hash_db ){
	def file = new File(hash_db)

	if (!file.exists()) {
	    throw new RuntimeException("File 'md5sum.txt' does not exist.")
	}
	if (!file.isFile()) {
	    throw new RuntimeException("Path 'md5sum.txt' is not a regular file.")
	}
	def hashMap = [:]
	file.eachLine { line ->
	    def parts = line.split(/\s+/).findAll { it.trim() != '' }
	    if (parts.size() >= 2) {
		def hash = parts[0]
		def filename = parts[1]
		hashMap[filename] = hash
	    }
	}
	return hashMap
}

def removeSuffix(String str, String suffix) {
    if (str == null || suffix == null) return str
        return str.endsWith(suffix) ? str[0..<str.length() - suffix.length()] : str
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
    //def hashes = '' //args.hash_db

    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    samplesheet_flavor = "modal"

    if (samplesheet_flavor == "flat"){
	    header = ['sample_id', 'normal_read1', 'normal_read2', 'tumor_read1', 'tumor_read2', 'normal_read1_hash', 'normal_read2_hash', 'tumor_read1_hash', 'tumor_read2_hash']
	}

    // modality is either FF, FFPE, BLOOD
    if (samplesheet_flavor == "modal"){
	    header = ['sample_id', 'normal_modality', 'tumor_modality']
	}

    csv_channel = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: header, skip: 1)


    if(samplesheet_flavor == 'flat'){
            // cast the types of the channel
            sample_pairs = csv_channel.map{
                row -> 
                [sample_id:row.sample_id, 
                normal_read1:file(row.normal_read1),
                normal_read2:file(row.normal_read2),
                tumor_read1:file(row.tumor_read1),
                tumor_read2:file(row.tumor_read2),
                normal_read1_hash:row.normal_read1_hash,
                normal_read2_hash:row.normal_read2_hash,
                tumor_read1_hash:row.tumor_read1_hash,
                tumor_read2_hash:row.tumor_read2_hash,
                ]
            }

            // explicitely stage inputs
            staged_sample_pairs = stage_fastq(sample_pairs)

            // still check the rows for the naming scheme, ignoring the hashes
	    def check_row = { row ->
		def sid = row[0]
		def expected_row = [
			"${sid}",
			"${sid}_N_1.fq.gz",
			"${sid}_N_2.fq.gz",
			"${sid}_T_1.fq.gz",
			"${sid}_T_2.fq.gz",
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
	    
	    // a flat channel of [sample_id, read] tuples, for maximum parallelism
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
                row -> [[row[0], row[1], row[2], "${row[0]}_N", row[5], row[6]],
                        [row[0], row[3], row[4], "${row[0]}_T", row[7], row[8]]
                        ]
                }
    } else if(samplesheet_flavor == 'modal'){

        csv_channel.subscribe{ row ->
          if (!(row.normal_modality in ["FFPE","BLOOD"])){
                throw new Exception("unknown modality ${row.tumor_modality} in row: ${row}")
          } else if (!(row.tumor_modality in ["FFPE","FF"])){
                throw new Exception("unknown modality ${row.tumor_modality} in row: ${row}")
          }
        }

        
        def input_dir_ = ''
        def input_dir = ''

	// Validate and normalize directory path
        if (!args.containsKey('input_dir')){
          input_dir_ = "${launchDir}"
        } else {
	  input_dir_ = args.input_dir ?: "${launchDir}"
	}

	if (input_dir_) {
	    def dirObj = new File(input_dir_)
	    if (!dirObj.exists()) {
		throw new IllegalArgumentException("Directory does not exist: ${dir}")
	    }
	    if (!dirObj.isDirectory()) {
		throw new IllegalArgumentException("Path is not a directory: ${dir}")
	    }
	    input_dir = dirObj.absolutePath
	}

	def hashes = args.hash_db


        assert hashes != null : "hashdb argument cant be null when running samplesheet_flavor: modal"
        def hashes_map = hash_db_to_dict(hashes)
	println "hashes map: ${hashes_map}"
        //throw new Exception("hashes map: ${hashes_map}")

        
        sample_reads_w_prefix = csv_channel.map{ row -> 
        [[ row.sample_id, 
                file("${input_dir}/${row.sample_id}_N_${row.normal_modality}_1.fq.gz"),
                file("${input_dir}/${row.sample_id}_N_${row.normal_modality}_2.fq.gz"),
                "${row.sample_id}_N",
                hashes_map["${row.sample_id}_N_${row.normal_modality}_1.fq.gz"],
                hashes_map["${row.sample_id}_N_${row.normal_modality}_2.fq.gz"],
                ],
        [ row.sample_id, 
                file("${input_dir}/${row.sample_id}_T_${row.tumor_modality}_1.fq.gz"),
                file("${input_dir}/${row.sample_id}_T_${row.tumor_modality}_2.fq.gz"),
                "${row.sample_id}_T",
                hashes_map["${row.sample_id}_T_${row.tumor_modality}_1.fq.gz"],
                hashes_map["${row.sample_id}_T_${row.tumor_modality}_2.fq.gz"],
                ]
        ]
        }.flatMap{ it -> [it[0], it[1]] }

        bam_pairings = csv_channel.map{ row -> 
                ["${row.sample_id}_N_${row.normal_modality}_1",
                "${row.sample_id}_T_${row.tumor_modality}_1"
                ]}

    } else {
        throw new Exception("unsupported samplesheet flavor: ${samplesheet_flavor}")
    }


    // output_file_prefix has to be calculated here, 
    // because it has to be known before starting the process as it it an output thereof
    fastp_out = fastp(sample_reads_w_prefix.unique())
    preprocessed_reads = fastp_out.preprocessed_reads
    integrity_check = fastp_out.integrity_check
    fastp_report_h = fastp_out.html
    fastp_report_j = fastp_out.json

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

    if(args.bwa_tool == 'bwa'){
    	bams = bwa_align(preprocessed_reads, args.refgenome, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt, bwa_idx.pac, bwa_idx.sa, args.cleanup_intermediate_files)
    } else if(args.bwa_tool == 'bwa2'){
        // bwa_idx.f0123 was has an f prefixed so, the file ending starts with a non-numeral
        bams = bwamem2_align(preprocessed_reads, args.refgenome, bwa_idx.f0123, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt_2bit_64, bwa_idx.pac, args.cleanup_intermediate_files)
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
    //bam_w_depth = bam_depth(bam_recalibrated, args.refgenome)
    bam_coverage = bam_coverage(bam_recalibrated)
    bams_w_stats = bam_stats(bam_recalibrated, args.refgenome)


    variant_out = variant_call(
            bam_pairings,
	    bam_recalibrated,
	    args.intervals,
	    args.panel_of_normals,
	    args.germline_resource,
	    args.refgenome,
	    args,
	    )

  emit:
    bam_pairs_w_idx = variant_out.bam_pairs_w_idx
    vcf = variant_out.vcf
    //bam_depth = bam_w_depth
    bam_coverage = bam_coverage
    bam_stats = bams_w_stats
    refgenome_index = refgenome_index
    refgenome_dict = refgenome_dict
    samplesheet = Channel.of(args.samplesheet)
    hash_db = Channel.of(args.hash_db)
    integrity_check = integrity_check
    fastp_report_h = fastp_report_h
    fastp_report_j = fastp_report_j
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}
