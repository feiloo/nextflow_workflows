process gatk_createsequencedictionary {
    // note, this process requires already sorted bams

    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    storeDir "$NEXTFLOW_STOREDIR"

    cpus { Math.max(1, Math.round(8 * (1 - ((1/4)*(task.attempt-1))))) }
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

    cpus { Math.max(1, Math.round(8 * (1 - ((1/4)*(task.attempt-1))))) }
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
    conda "bioconda::gatk4=4.4.0.0 bioconda::samtools=1.17"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0 quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    cpus { Math.max(1, Math.round(26 * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory {(80+(80 * (task.attempt-1))).GB}

    input:
        path(bam)
        path(intervals)
        val(library_type)

    output:
      path("out/${bam}"), emit: marked_bams
      //path("${bam.getSimpleName()}_metrics.txt"), emit: metrics

    script:
    // write to output directory to avoid filename collision but keeping the filename in the channels the same
    //--METRICS_FILE ${bam.getSimpleName()}_metrics.txt
    """
    mkdir out
    mkdir tmp

    if [[ "${library_type}" == "wes" ]]
    then
        genomic_region="--intervals ${intervals}"
        samtools index "${bam}" -@ ${task.cpus} -o "${bam}.bai"
    elif [[ "${library_type}" == "wgs" ]]
    then
        genomic_region=""
    fi

    gatk MarkDuplicatesSpark \\
        --input ${bam} \\
        --java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G" \\
	--conf 'spark.local.dir=tmp' \\
	--conf 'spark.executor.cores=${task.cpus}' \\
	--conf 'spark.executor.memory=24g' \\
	--conf 'spark.driver.memory=2g' \\
	--tmp-dir tmp \\
    	--spark-master local[${task.cpus}] \\
	--output out/${bam} \\
        \$genomic_region 
        

    """
}

process gatk_set_tags {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'


    cpus { Math.max(1, Math.round(8 * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory {Math.min(80, 80+(80 * (task.attempt-1))).GB}

    input:
        path(bam)
	path(refgenome)
	val(cleanup_intermediate_files)

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

    if [[ "${cleanup_intermediate_files}" == 'true' ]]; then
      rm ${bam}
    fi
    """
}

process gatk_baserecalibrator {
    conda "bioconda::gatk4=4.4.0.0 bioconda::samtools=1.17"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0 quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    cpus { Math.max(1, Math.round(12 * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4


    memory {Math.min(96, 20+(14 * (task.attempt-1))).GB}


    input:
        path(bamfile)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)
	path(known_sites)
	path(known_sites_index)
        path(intervals)
        val(library_type)

    output:
        path("${bamfile.getSimpleName()}_recal_data.table")

    script:
    """
    mkdir tmp

    if [[ "${library_type}" == "wes" ]]
    then
        genomic_region="--intervals ${intervals}"
        samtools index "${bamfile}" -@ ${task.cpus} -o "${bamfile}.bai"
    elif [[ "${library_type}" == "wgs" ]]
    then
        genomic_region=""
    fi

    gatk BaseRecalibrator \\
	--java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G -XX:ParallelGCThreads=2" \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--known-sites ${known_sites} \\
	--output ${bamfile.getSimpleName()}_recal_data.table \\
        \$genomic_region
    """

}

process gatk_apply_bqsr {
    conda "bioconda::gatk4=4.4.0.0 bioconda::samtools=1.17"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0 quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    memory {Math.min(96, 20+(14 * (task.attempt-1))).GB}

    cpus { Math.max(1, Math.round(12 * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4


    input:
        tuple path(bamfile), path(bam_recal_data)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)
        path(intervals)
        val(library_type)

    output:
        path("out/${bamfile.getSimpleName()}.bam")

    script:
    """
    mkdir out
    mkdir tmp

    if [[ "${library_type}" == "wes" ]]
    then
        genomic_region="--intervals ${intervals}"
        samtools index "${bamfile}" -@ ${task.cpus} -o "${bamfile}.bai"
    elif [[ "${library_type}" == "wgs" ]]
    then
        genomic_region=""
    fi

    gatk ApplyBQSR \\
        --java-options "-Djava.io.tmpdir=tmp -Xms24G -Xmx24G -XX:ParallelGCThreads=2" \\
        --input ${bamfile} \\
	--reference ${refgenome} \\
	--bqsr-recal-file ${bam_recal_data} \\
	--output out/${bamfile.getSimpleName()}.bam \\
        \$genomic_region
    """

}

// test
workflow {
  filec = Channel.fromPath(params.bam)
  refg = Channel.fromPath('/data/reference/clc_refseqs/Homo_sapiens_sequence_hg38_no_alt_analysis_set.fa')
  gatk_markduplicates(filec, refg)
}
