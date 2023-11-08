process gatk_markduplicates {
    // note, this process requires already sorted bams

    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        tuple val(sample_id), path(bam)
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

// test
workflow {
  filec = Channel.fromPath(params.bam)
  refg = Channel.fromPath('/data/reference/clc_refseqs/Homo_sapiens_sequence_hg38_no_alt_analysis_set.fa')
  fc2 = filec.map{it->['sample_s1', it]}
  gatk_markduplicates(fc2, refg)
}
