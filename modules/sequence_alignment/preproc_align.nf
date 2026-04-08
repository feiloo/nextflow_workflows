process fastp {
    conda "bioconda::fastp=1.0.1 bioconda::bcftools=1.17"
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'

    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    cache 'lenient'

    //cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    cpus 12
    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4
    // highest compression by default
    ext compression_level: "9"

    input:
    tuple val(sample_id), path(read1), path(read2), val(output_file_prefix), val(expected_checksum1), val(expected_checksum2)

    output:
    tuple path("out/${read1.getSimpleName()}.fq.gz"), path("out/${read2.getSimpleName()}.fq.gz"), emit: preprocessed_reads
    path("${output_file_prefix}_fastp.html"), emit: html
    path("${output_file_prefix}_fastp.json"), emit: json
    path("${output_file_prefix}_check_process_results.txt"), emit: integrity_check

    script:


    """
    set -euo pipefail

    GZ_COMPRESSION_LEVEL=${task.ext.compression_level ?: "9"} # out of 1 to 9

    mkdir -p out
    CHECKSUM_TYPE="md5sum"

    fastp_preprocessing() {
      fastp \\
	--in1 ${read1} \\
	--in2 ${read2} \\
	--out1  out/${read1.getSimpleName()}.fq.gz \\
	--out2  out/${read2.getSimpleName()}.fq.gz \\
	-z \$GZ_COMPRESSION_LEVEL \\
	--thread ${task.cpus} \\
	--json ${output_file_prefix}_fastp.json \\
	--html ${output_file_prefix}_fastp.html \\
	2> ${output_file_prefix}.fastp.log
	}


    check_integrity() {
	    NAME=\$1
	    NAMEOUT="\$NAME"_out
	    EXPECTED_CHECKSUM=\$2

	    mkdir "\$NAME"_out
            \$CHECKSUM_TYPE \$NAME | awk '{print \$1}' > \$NAMEOUT/"\$CHECKSUM_TYPE"_result &
	    CHECKSUM_PID=\$!

	    gzip -t --keep \$NAME > \$NAMEOUT/gzip_result &
	    GZIP_PID=\$!
	    bgzip -t -@ 8 --keep \$NAME > \$NAMEOUT/bgzip_result
	    BGZIP_PID=\$!

	    wait \$CHECKSUM_PID

	    CHECKSUM_RESULT=\$(cat \$NAMEOUT/"\$CHECKSUM_TYPE"_result)
	    if [ "\$CHECKSUM_RESULT" != "\$EXPECTED_CHECKSUM" ]; then
		echo "error, \$CHECKSUM_TYPE of \$NAME doesnt match: \$CHECKSUM_RESULT  \$EXPECTED_CHECKSUM"
		exit 250

	    fi 
	    #MD5SUM_PID=\$!
	    wait \$GZIP_PID

	    if grep -qw 'trailing' \$NAMEOUT/gzip_result && grep -qw 'garbage' \$NAMEOUT/gzip_result; then
		echo "error, found trailing garbage in \$NAME according to gzip. the file is possibly truncated"
		exit 251
	    fi

	    # todo add verficiation of bgzip too here
	    wait \$BGZIP_PID

	    echo -e "check process results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/"\$CHECKSUM_TYPE"_result >> \$NAMEOUT/check_process_results
	    echo -e "\ngzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/gzip_result >> \$NAMEOUT/check_process_results
	    echo -e "\nbgzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/bgzip_result >> \$NAMEOUT/check_process_results

    }

    fastp_preprocessing &
    FASTP_PID=\$!

    check_integrity "${read1}" "${expected_checksum1}" &
    CHECK1=\$!
    check_integrity "${read2}" "${expected_checksum2}" &
    CHECK2=\$!

    wait \$CHECK1
    wait \$CHECK2
    wait \$FASTP_PID

    cat ${read1}_out/check_process_results ${read2}_out/check_process_results > ${output_file_prefix}_check_process_results.txt

    """
}

process fastp_bwa {
    conda "bioconda::fastp=1.0.1 bioconda::bcftools=1.17 bioconda::bwa-mem2=2.2.1 samtools=1.16.1"
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'
    //container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0'

    memory "${144 * (1+0.3*task.attempt)} GB"
    cache 'lenient'

    cpus 26 // + 12

    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4
    // highest compression by default
    ext compression_level: "9"

    input:
    tuple val(sample_id), path(read1), path(read2), val(output_file_prefix), val(expected_checksum1), val(expected_checksum2)
    path(refgenome)
    path("${refgenome}.0123")
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt.2bit.64")
    path("${refgenome}.pac")


    output:
    path("${read1.getSimpleName()}.bam"), emit: bam
    path("${output_file_prefix}_fastp.html"), emit: html
    path("${output_file_prefix}_fastp.json"), emit: json
    path("${output_file_prefix}_check_process_results.txt"), emit: integrity_check

    script:

    // see https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-dictionary/Read_groups.md for more info
    // example: def read_group_info = "@RG\tID:SAMN10471711\tSM:SAMN10471711\tLB:SAMN10471711\tPL:ILLUMINA"
    // the real read group data has to come from the illumina casava 1.8 format lines from the fastq

    // alternative, use picard AddOrReplaceReadGroups
    // def read_group_info = "@RG\tID:SAMN10471711\tSM:SAMN10471711\tLB:SAMN10471711\tPL:ILLUMINA"
    // format: "instrument_name":"run_id":"flowcell_id":"flowcell_lane":"tile number within flowcell":"x coordinate of cluster within tile":"y coordinate of cluster within tile"<space>"member of a pair 1 or 2, paired end or mate pair reads only":"Y if read is filtered did not pass, N otherwise":"control bits":"index sequence"
    // example 
    //@EAS139:136:FC706VJ:1:2104:15343:197393 1:N:0:ATCACG

    // todo, write a proper readgroup parser/detector

    // workaround this by
    // assuming the inputs are already demultiplexed, so just add a stub readgroup
    def read_group_identifier = "RunId1"

    // this is just fake and testdata for now
    def flowcell_barcode = "FC706VJ"
    def flowcell_lane = "1"
    def sample_barcode = "${read1.getSimpleName()}"

    def platform_unit = "${flowcell_barcode}.${flowcell_lane}.${sample_barcode}"
    def sample_name = "${read1.getSimpleName()}"
    
    // one of like: ILLUMINA, SOLID, LS454, HELICOS, PACBIO
    def platform_technology = "ILLUMINA"

    def library_prep_identifier = "prep1"

    def read_group_info = "@RG\\tID:${read_group_identifier}\\tPL:${platform_technology}\\tLB:${library_prep_identifier}\\tPU:${platform_unit}\\tSM:${sample_name}"

    // todo: remove trailing read-numbers like "T_1.bam" from bam since bams combine both reads here


    """
    set -euo pipefail

    mkdir -p out
    CHECKSUM_TYPE="md5sum"

    fastp_preprocessing() {
      fastp \\
	--in1 ${read1} \\
	--in2 ${read2} \\
	--stdout \\
	--thread 4 \\
	--dedup \\
	--dup_calc_accuracy 6 \\
	--json ${output_file_prefix}_fastp.json \\
	--html ${output_file_prefix}_fastp.html \\
	2> ${output_file_prefix}.fastp.log \\
	| bwa-mem2 mem -R "${read_group_info}" -t 20 ${refgenome} - \\
	| samtools view --bam --threads 2 -o ${read1.getSimpleName()}.bam -
	}


    check_integrity() {
	    NAME=\$1
	    NAMEOUT="\$NAME"_out
	    EXPECTED_CHECKSUM=\$2

	    mkdir "\$NAME"_out
            \$CHECKSUM_TYPE \$NAME | awk '{print \$1}' > \$NAMEOUT/"\$CHECKSUM_TYPE"_result &
	    CHECKSUM_PID=\$!

	    gzip -t --keep \$NAME > \$NAMEOUT/gzip_result &
	    GZIP_PID=\$!
	    bgzip -t -@ 8 --keep \$NAME > \$NAMEOUT/bgzip_result
	    BGZIP_PID=\$!

	    wait \$CHECKSUM_PID

	    CHECKSUM_RESULT=\$(cat \$NAMEOUT/"\$CHECKSUM_TYPE"_result)
	    if [ "\$CHECKSUM_RESULT" != "\$EXPECTED_CHECKSUM" ]; then
		echo "error, \$CHECKSUM_TYPE of \$NAME doesnt match: \$CHECKSUM_RESULT  \$EXPECTED_CHECKSUM"
		exit 250

	    fi 
	    #MD5SUM_PID=\$!
	    wait \$GZIP_PID

	    if grep -qw 'trailing' \$NAMEOUT/gzip_result && grep -qw 'garbage' \$NAMEOUT/gzip_result; then
		echo "error, found trailing garbage in \$NAME according to gzip. the file is possibly truncated"
		exit 251
	    fi

	    # todo add verficiation of bgzip too here
	    wait \$BGZIP_PID

	    echo -e "check process results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/"\$CHECKSUM_TYPE"_result >> \$NAMEOUT/check_process_results
	    echo -e "\ngzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/gzip_result >> \$NAMEOUT/check_process_results
	    echo -e "\nbgzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/bgzip_result >> \$NAMEOUT/check_process_results

    }

    fastp_preprocessing &
    FASTP_PID=\$!

    check_integrity "${read1}" "${expected_checksum1}" &
    CHECK1=\$!
    check_integrity "${read2}" "${expected_checksum2}" &
    CHECK2=\$!

    wait \$CHECK1
    wait \$CHECK2
    wait \$FASTP_PID

    cat ${read1}_out/check_process_results ${read2}_out/check_process_results > ${output_file_prefix}_check_process_results.txt

    """
}

process bwamem2_index_refgenome {
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    memory '180 GB'
    storeDir "$NEXTFLOW_STOREDIR/bwamem2_indices"

    input:
    path(refgenome)

    output:

    // file endings are optionally prefix with an f so they start with an non-numerals
    // as a convention
    path("${refgenome}.0123"), emit: f0123
    path("${refgenome}.amb"), emit: amb
    path("${refgenome}.ann"), emit: ann
    path("${refgenome}.bwt.2bit.64"), emit: bwt_2bit_64
    path("${refgenome}.pac"), emit: pac

    script:
    """
    bwa-mem2 index ${refgenome}
    """
}

process bwamem2_align {
    conda "bioconda::bwa-mem2=2.2.1 samtools=1.16.1"
    container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0'

    memory "${128 * (1+0.3*task.attempt)} GB"
    cpus 26

    input:
    tuple path(read1), path(read2)
    path(refgenome)
    path("${refgenome}.0123")
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt.2bit.64")
    path("${refgenome}.pac")
    val(cleanup_intermediate_files)

    output:
    path("${read1.getSimpleName()}.bam")

    script:
    // see https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-dictionary/Read_groups.md for more info
    // example: def read_group_info = "@RG\tID:SAMN10471711\tSM:SAMN10471711\tLB:SAMN10471711\tPL:ILLUMINA"
    // the real read group data has to come from the illumina casava 1.8 format lines from the fastq

    // alternative, use picard AddOrReplaceReadGroups
    // def read_group_info = "@RG\tID:SAMN10471711\tSM:SAMN10471711\tLB:SAMN10471711\tPL:ILLUMINA"
    // format: "instrument_name":"run_id":"flowcell_id":"flowcell_lane":"tile number within flowcell":"x coordinate of cluster within tile":"y coordinate of cluster within tile"<space>"member of a pair 1 or 2, paired end or mate pair reads only":"Y if read is filtered did not pass, N otherwise":"control bits":"index sequence"
    // example 
    //@EAS139:136:FC706VJ:1:2104:15343:197393 1:N:0:ATCACG

    // todo, write a proper readgroup parser/detector

    // workaround this by
    // assuming the inputs are already demultiplexed, so just add a stub readgroup
    def read_group_identifier = "RunId1"

    // this is just fake and testdata for now
    def flowcell_barcode = "FC706VJ"
    def flowcell_lane = "1"
    def sample_barcode = "${read1.getSimpleName()}"

    def platform_unit = "${flowcell_barcode}.${flowcell_lane}.${sample_barcode}"
    def sample_name = "${read1.getSimpleName()}"
    
    // one of like: ILLUMINA, SOLID, LS454, HELICOS, PACBIO
    def platform_technology = "ILLUMINA"

    def library_prep_identifier = "prep1"

    def read_group_info = "@RG\\tID:${read_group_identifier}\\tPL:${platform_technology}\\tLB:${library_prep_identifier}\\tPU:${platform_unit}\\tSM:${sample_name}"

    // todo: remove trailing read-numbers like "T_1.bam" from bam since bams combine both reads here

    """
    bwa-mem2 mem -R "${read_group_info}" -t ${task.cpus} ${refgenome} ${read1} ${read2} | samtools view --bam --threads ${task.cpus} -o ${read1.getSimpleName()}.bam
    """
}

workflow bwamem2_tumor_only {
  take:
    args
  main:
    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)

    i1 = Channel.fromPath("${args.refgenome}.0123")
    i2 = Channel.fromPath("${args.refgenome}.bwt.2bit.64")
    i3 = Channel.fromPath("${args.refgenome}.alt")
    i4 = Channel.fromPath("${args.refgenome}.ann")
    i5 = Channel.fromPath("${args.refgenome}.pac")

    bwamem2_align(csv_channel, args.refgenome, i1, i2, i3, i4, i5)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  bwamem2_tumor_only(args)

}



