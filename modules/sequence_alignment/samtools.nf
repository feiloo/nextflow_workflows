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
    cpus 16

    input:
    path(samfile)
    val(cleanup_intermediate_files)

    output:
    path("${samfile.getBaseName()}.bam")

    script:
    """
    samtools view "${samfile}" --bam --threads ${task.cpus} -o "${samfile.getBaseName()}.bam"
    if [[ "${cleanup_intermediate_files}" == 'true' ]]; then
          rm ${samfile}
    fi
    """
}

process sort_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'


    //samtools sort uses by default 768MB per thread
    // so we set a limit according to it
    cpus 16
    memory "20 GB"

    input:
    path(bamfile)

    output:
    path("out/${bamfile}")

    script:
    """
    mkdir out
    samtools sort "${bamfile}" -@ ${task.cpus} -o "out/${bamfile}"
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
    """
    samtools index "${bamfile}" -@ ${task.cpus} -o "${bamfile}.bai"
    """
}

process bam_stats {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 4

    memory "1 GB"

    input:
    path(bamfile)
    path(refgenome)

    output:
    path("${bamfile.getSimpleName()}_samstats")

    script:
    n_cpus = 4
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

    memory "8 GB"
    input:
    path(bamfile)
    path(refgenome)

    output:
    path("${bamfile.getSimpleName()}_samdepth")

    script:
    """
    samtools depth -@ ${task.cpus} "${bamfile}" > "${bamfile.getSimpleName()}_samdepth"
    """
}

process bam_coverage {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus 1

    memory "1 GB"
    input:
    path(bamfile)

    output:
    path("${bamfile.getSimpleName()}_samcov")

    script:
    """
    samtools coverage "${bamfile}" -o "${bamfile.getSimpleName()}_samcov"
    """
}
