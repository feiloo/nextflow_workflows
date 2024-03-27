process bwa_index_refgenome {
    conda "bioconda::bwa=0.7.17"
    //container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    storeDir "$NEXTFLOW_STOREDIR"

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    memory '120 GB'

    input:
    path(refgenome)

    output:

    // file endings are optionally prefix with an f so they start with an non-numerals
    // as a convention
    path("${refgenome}.amb"), emit: amb
    path("${refgenome}.ann"), emit: ann
    path("${refgenome}.bwt"), emit: bwt
    path("${refgenome}.pac"), emit: pac
    path("${refgenome}.sa"), emit: sa

    script:
    """
    bwa index ${refgenome}
    """
}

process bwa_align {
    conda "bioconda::bwa=0.7.17 samtools=1.16.1"
    //container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    publishDir '/PAT-Sequenzer/NEB_FFPE_WGS_30-01-2024/nextflow_outputs/other_bams', mode: 'copy', overwrite: true

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple path(read1), path(read2)
    path(refgenome)
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt")
    path("${refgenome}.pac")
    path("${refgenome}.sa")
    val(cleanup_intermediate_files)


    output:
    path("${read1.getSimpleName()}.bam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    // todo: fix read groups:
    // def read_group_info = "read group"
    // bwa mem -R ${read_group_info} -t $n_cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam


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

    // bwa mem -t $n_cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam
    """
    bwa mem -R "${read_group_info}" -t ${task.cpus} ${refgenome} ${read1} ${read2} | samtools view --bam --threads ${task.cpus} -o ${read1.getSimpleName()}.bam
    if [[ "${args.cleanup_intermediate_files}" == 'true' ]]; then
      rm ${read1} && rm ${read2}
    fi
    """
}


