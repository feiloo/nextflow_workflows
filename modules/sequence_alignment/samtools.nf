process index_fasta {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    storeDir "$NEXTFLOW_STOREDIR"

    input:
    path(refgenome)

    output:
    path("${refgenome}.fai"), emit: fasta_index

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    samtools faidx ${refgenome} --output ${refgenome}.fai
    """
}

process sam_to_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    memory "10 GB"

    input:
    path(samfile)
    val(cleanup_intermediate_files)

    output:
    path("${samfile.getBaseName()}.bam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    samtools view "${samfile}" --bam --threads $n_cpus -o "${samfile.getBaseName()}.bam"
    if [[ "${cleanup_intermediate_files}" == 'true' ]]; then
          rm ${samfile}
    fi
    """
}

process sort_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path(bamfile)

    output:
    path("${bamfile}.bam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    mv ${bamfile} unsorted_${bamfile}
    samtools sort "unsorted_${bamfile}" -@ $n_cpus -o "${bamfile}.bam"
    """
}

process index_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 8

    input:
    path(bamfile)

    output:
    path("${bamfile}.bai")

    script:
    n_cpus = 8
    """
    samtools index "${bamfile}" -@ $n_cpus -o "${bamfile}.bai"
    """
}

process bam_stats {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 8

    memory "1 GB"

    input:
    path(bamfile)
    path(refgenome)

    output:
    path("${bamfile.getSimpleName()}_samstats")

    script:
    n_cpus = 8
    """
    samtools stats "${bamfile}" \\
    	--reference ${refgenome} \\
	-@ $n_cpus \\
	> "${bamfile.getSimpleName()}_samstats"
    """
}

process bam_depth {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 8

    memory "1 GB"
    input:
    path(bamfile)
    path(refgenome)

    output:
    path("${bamfile.getSimpleName()}_samdepth")

    script:
    n_cpus = 8
    """
    samtools depth "${bamfile}" > "${bamfile.getSimpleName()}_samdepth"
    """
}

process bam_coverage {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 8

    memory "1 GB"
    input:
    path(bamfile)

    output:
    path("${bamfile.getSimpleName()}_samcov")

    script:
    n_cpus = 8
    """
    samtools coverage "${bamfile}" -o "${bamfile.getSimpleName()}_samcov"
    """
}
