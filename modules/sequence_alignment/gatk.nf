process gatk_createsequencedictionary {
    // note, this process requires already sorted bams

    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    storeDir "$NEXTFLOW_STOREDIR"

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory "56 GB"

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

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

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

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory {Math.min(56, 30+(14 * (task.attempt-1))).GB}

    input:
        path(bam)

    output:
      path("out/${bam}"), emit: marked_bams
      //path("${bam.getSimpleName()}_metrics.txt"), emit: metrics

    script:
    // write to output directory to avoid filename collision but keeping the filename in the channels the same
    //--METRICS_FILE ${bam.getSimpleName()}_metrics.txt

    """
    mkdir out
    mkdir tmp
    gatk MarkDuplicatesSpark \\
        --input ${bam} \\
	--java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G" \\
	--conf 'spark.local.dir=tmp' \\
	--conf 'spark.executor.cores=${task.cpus}' \\
	--conf 'spark.executor.memory=24g' \\
	--conf 'spark.driver.memory=2g' \\
	--tmp-dir tmp \\
    	--spark-master local[${task.cpus}] \\
	--output out/${bam}

    """
}

process gatk_set_tags {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'


    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory {Math.min(56, 35+(14 * (task.attempt-1))).GB}

    input:
        path(bam)
	path(refgenome)

    output:
      path("out/${bam}"), emit: tagged_bams

    script:
    """
    mkdir out
    mkdir tmp
    gatk SetNmMdAndUqTags \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
    	--INPUT ${bam} \\
	--REFERENCE_SEQUENCE ${refgenome} \\
	--OUTPUT out/${bam}
    """
}

process gatk_baserecalibrator {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4


    memory {Math.min(56, 20+(14 * (task.attempt-1))).GB}


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
    mkdir tmp
    gatk BaseRecalibrator \\
	--java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G -XX:ParallelGCThreads=2" \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--known-sites ${known_sites} \\
	--output ${bamfile.getSimpleName()}_recal_data.table
    """

}

process gatk_apply_bqsr {
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    memory {Math.min(56, 20+(14 * (task.attempt-1))).GB}

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4


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
    mkdir tmp
    gatk ApplyBQSR \\
        --java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G -XX:ParallelGCThreads=2" \\
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
