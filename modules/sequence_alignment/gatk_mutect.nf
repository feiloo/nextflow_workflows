// implement best-practice somatic variant calling
// reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2


process gatk_mutect {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(sample_bam)
	path(intervals)
	path(panel_of_normals)
	path(germline_resource)
	path(refgenome)

    output:
      path("${sample_bam.getSimpleName()}_unfiltered.vcf"), emit: vcf
      path("${sample_bam.getSimpleName()}_f1r2_data.tar.gz"), emit: f1r2_data

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    gatk Mutect2 \\
    	--native-pair-hmm-threads ${n_cpus} \\
	--INPUT ${sample_bam} \\
	-L ${intervals} \\
	-pon ${panel_of_normals} \\
	-germline-resource ${germline_resource} \\
        --REFERENCE ${refgenome} \\
	--f1r2-tar-gz "${sample_bam.getSimpleName()}_f1r2_data.tar.gz" \\
	--OUTPUT ${sample_bam.getSimpleName()}_unfiltered.vcf}
    """
}

process gatk_learn_readorientationmodel {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(f1r2_data)

    output:
      path("${f1r2_data.getSimpleName()}_model.tar.gz"), emit: orientation_model

    script:
    """
    gatk LearnReadOriantationModel \\
	--INPUT ${f1r2_data} \\
	--OUTPUT "${f1r2_data.getSimpleName()}_model.tar.gz"
    """
}


process gatk_getpileupsummaries {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(sample_bam)
        path(sample_vcf)
        path(genomic_intervals)
	// prepared from like gnomAD with pop-freqs in the info field
        path(variant_frequency_vcf)

    output:
      path("${sample_vcf.getSimpleName()}_pileup.table"), emit: vcf

    script:
    """
    gatk GetPileupSummaries \\
	--INPUT ${f1r2_data} \\
	-L ${genomic_intervals} \\
	-V ${variant_frequency_vcf} \\
	--OUTPUT "${sample_vcf.getSimpleName()}_pileup.table"
    """
}

process gatk_calculate_contamination {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        tuple path(tumor_pileups), path(normal_pileups)
        path(segments_table)

    output:
      path("${tumor_pileups.getSimpleName()}_contamination.table"), emit: contamination_table

    script:
    """
    gatk CalculateContamination \\
	--INPUT ${tumor_pileups} \\
	-matched ${normal_pileups} \\
	--OUTPUT "${tumor_pileups.getSimpleName()}_contamination.table"
    """
}

process gatk_filter_calls {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(sample_vcf)
        path(segments_table)
        path(contamination_table)
	path(orientation_model)

    output:
      path("${sample_vcf.getSimpleName()}_filtered.vcf"), emit: vcf

    script:
    """
    gatk FilterMutectCalls \\
	--variant ${sample_vcf} \\
	--tumor-segmentation ${segments_table} \\
	--ob-priors ${orientation_model} \\
	--OUTPUT "${sample_vcf.getSimpleName()}_filtered.vcf"
    """
}


workflow sequence_alignment {
  take:
    args
  main:
    mut = gatk_mutect(sample_bam, intervals, panel_of_normals, germline_resourse, args.refgenome)
    om = gatk_learn_readorientationmodel(mut.f1r2_data).orientation_model

    # for tumors and for normals
    pileups = gatk_getpileupsummaries(sample_bam, mut.vcf, genomic_intervals, variant_frequency_vcf)

    matched_pileups = tumor_pileups.cross(normal_pileups)
    contamination = gatk_calculate_contamination(matched_pileups, segments_table)
    filtered_vcf = gatk_filter_calls(mut.vcf, segments_table, contamination, om)

  emit:
    vcf = filtered_vcf
}
