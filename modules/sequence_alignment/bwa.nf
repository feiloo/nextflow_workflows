process bwa_index_refgenome {
    conda "bioconda::bwa=0.7.17"
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    storeDir "$NEXTFLOW_STOREDIR"

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    //memory '120 GB'

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
    conda "bioconda::bwa=0.7.17"
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    input:
    tuple path(read1), path(read2)
    path(refgenome)
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt")
    path("${refgenome}.pac")
    path("${refgenome}.sa")


    output:
    path("${read1.getSimpleName()}.sam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    bwa mem -t $n_cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam
    """
}


