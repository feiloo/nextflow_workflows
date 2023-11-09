process gatk_markduplicates {
    // note, this process requires already sorted bams

    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(bam)
	path(refgenome)

    output:
      path("out/${bam}")
      path("${bam.getSimpleName()}_metrics.txt")

    script:
    // write to output directory to avoid filename collision but keeping the filename in the channels the same
    """
    mkdir out
    gatk MarkDuplicates \\
        --INPUT ${bam} \\
	--OUTPUT out/${bam} \\
	--METRICS_FILE ${bam.getSimpleName()}_metrics.txt
    """
}

process gatk_baserecalibrator {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    input:
        path(bamfile)
	path(refgenome)
	path(known_sites)

    output:
        path("${bamfile.getSimpleName()}_recal_data.table")

    script:
    """
    gatk BaseRecalibrator \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--output ${bamfile.getSimpleName()}_recal_data.table
    """

}

// test
workflow {
  filec = Channel.fromPath(params.bam)
  refg = Channel.fromPath('/data/reference/clc_refseqs/Homo_sapiens_sequence_hg38_no_alt_analysis_set.fa')
  gatk_markduplicates(filec, refg)
}
