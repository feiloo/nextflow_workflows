process gatk_createsequencedictionary {
    // note, this process requires already sorted bams

    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    storeDir "$NEXTFLOW_STOREDIR"

    input:
	path(refgenome)

    output:
      path("${refgenome.getSimpleName()}.dict"), emit: refgenome_dict

    script:
    // write to output directory to avoid filename collision but keeping the filename in the channels the same
    """
    mkdir out
    gatk CreateSequenceDictionary \\
        --REFERENCE ${refgenome} \\
	--OUTPUT ${refgenome.getSimpleName()}.dict
    """
}

process gatk_indexfeaturefile {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    storeDir "$NEXTFLOW_STOREDIR"

    input:
        path(known_sites)

    output:
        path("${known_sites}.tbi"), emit: known_sites_index

    script:
    """
    gatk IndexFeatureFile \\
        --input ${known_sites} \\
	--output ${known_sites}.tbi
    """

}


process gatk_markduplicates {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory "60 GB"

    input:
        path(bam)

    output:
      path("out/${bam}"), emit: marked_bams
      //path("${bam.getSimpleName()}_metrics.txt"), emit: metrics

    script:
    // write to output directory to avoid filename collision but keeping the filename in the channels the same
    n_cpus = Runtime.runtime.availableProcessors()
    //--METRICS_FILE ${bam.getSimpleName()}_metrics.txt

    """
    mkdir out
    gatk MarkDuplicatesSpark \\
        --input ${bam} \\
    	--spark-master local[$n_cpus] \\
	--output out/${bam}

    """
}

process gatk_set_tags {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory "46 GB"

    input:
        path(bam)
	path(refgenome)
	val(cleanup_intermediate_files)

    output:
      path("out/${bam}"), emit: tagged_bams

    script:
    """
    mkdir out
    gatk SetNmMdAndUqTags \\
    	--INPUT ${bam} \\
	--REFERENCE_SEQUENCE ${refgenome} \\
	--OUTPUT out/${bam}

    if [[ "${cleanup_intermediate_files}" == 'true' ]]; then
      rm ${bam}
    fi
    """
}

process gatk_baserecalibrator {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    memory "46 GB"

    input:
        path(bamfile)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)
	path(known_sites)
	path(known_sites_index)

    output:
        path("${bamfile.getSimpleName()}_recal_data.table")

    script:
    """
    gatk BaseRecalibrator \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--known-sites ${known_sites} \\
	--output ${bamfile.getSimpleName()}_recal_data.table
    """

}

process gatk_apply_bqsr {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    memory "46 GB"

    input:
        tuple path(bamfile), path(bam_recal_data)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
        path("out/${bamfile.getSimpleName()}.bam")

    script:
    """
    mkdir out
    gatk ApplyBQSR \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--bqsr-recal-file ${bam_recal_data} \\
	--output out/${bamfile.getSimpleName()}.bam
    """

}

// test
workflow {
  filec = Channel.fromPath(params.bam)
  refg = Channel.fromPath('/data/reference/clc_refseqs/Homo_sapiens_sequence_hg38_no_alt_analysis_set.fa')
  gatk_markduplicates(filec, refg)
}
