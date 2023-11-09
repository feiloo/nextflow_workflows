process sam_to_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path(samfile)

    output:
    path("${samfile.getBaseName()}.bam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    samtools view ${samfile} --bam --threads $n_cpus -o ${samfile.getBaseName()}.bam
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

    input:
    path(bamfile)

    output:
    path("${bamfile}.bai")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    samtools index "${bamfile}" -@ $n_cpus -o "${bamfile}.bai"
    """
}

process bam_stats {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path(bamfile)
    path(refgenome)

    output:
    path("${bamfile}.bai")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    samtools stats "${bamfile}" \\
    	--reference ${refgenome} \\
	-@ $n_cpus
    """
}
