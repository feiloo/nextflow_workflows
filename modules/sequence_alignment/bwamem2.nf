process bwamem2_index_refgenome {
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    memory '120 GB'
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
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    memory "250 GB"

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
    path("${read1.getSimpleName()}.sam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
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

    """
    bwa-mem2 mem -R "${read_group_info}" -t $n_cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam

    if [[ "${cleanup_intermediate_files}" == 'true' ]]; then
      rm ${read1} && rm ${read2}
    fi
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

