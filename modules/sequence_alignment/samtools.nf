process sam_to_bam {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path(samfile)

    output:
    path("${samfile.getBaseName()}.bam")

    script:
    """
    samtools view ${samfile} --bam --threads 12 -o ${samfile.getBaseName()}.bam
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
    """
    mv ${bamfile} unsorted_${bamfile}
    samtools sort "unsorted_${bamfile}" -@ 12 -o "${bamfile}.bam"
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
    """
    samtools index "${bamfile}" -@ 12 -o "${bamfile}.bai"
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
    """
    samtools stats "${bamfile}" \\
    	--reference ${refgenome} \\
	-@ 12
    """
}
