// implement best-practice somatic variant calling
// reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2

include { index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

process gatk_mutect {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        tuple path(tumor_bam), path(normal_bam)
	path(intervals)
	path(panel_of_normals)
	path(panel_of_normals_index)
	path(germline_resource)
	path(germline_resource_index)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      path("${tumor_bam.getSimpleName()}_unfiltered.vcf"), emit: vcf
      path("${tumor_bam.getSimpleName()}_f1r2_data.tar.gz"), emit: f1r2_data

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    def normal_samplename = "${normal_bam.getSimpleName()}"
    """
    gatk Mutect2 \\
    	--native-pair-hmm-threads ${n_cpus} \\
	--input ${tumor_bam} \\
	--input ${normal_bam} \\
	-normal ${normal_bam.getSimpleName()} \\
	-L ${intervals} \\
	-pon ${panel_of_normals} \\
	-germline-resource ${germline_resource} \\
        --reference ${refgenome} \\
	--f1r2-tar-gz "${tumor_bam.getSimpleName()}_f1r2_data.tar.gz" \\
	--output "${tumor_bam.getSimpleName()}_unfiltered.vcf"
    """
}

process gatk_getsamplename {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
        path(sample_bam)

    output:
      val(header_samplename)

    script:
    """
    gatk GetSampleName \\
	--INPUT ${sample_bam}
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
    gatk LearnReadOrientationModel \\
	--input "${f1r2_data}" \\
	--output "${f1r2_data.getSimpleName()}_model.tar.gz"
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
        path(variant_frequency_vcf_index)

    output:
      path("${sample_vcf.getSimpleName()}_pileup.table"), emit: vcf

    script:
    """
    gatk GetPileupSummaries \\
	-I ${sample_bam} \\
	-L ${genomic_intervals} \\
	-V ${variant_frequency_vcf} \\
	-O "${sample_vcf.getSimpleName()}_pileup.table"
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

process create_pon_db {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
    	path(pon_vcfs)
    	path(refgenome)
    	path(intervals)

    output:
        path("pon_genomics_db")

    script:
    """
    gatk GenomicsDBImport \\
        -R ${refgenome} \\
        --genomicsdb-workspace-path pon_genomics_db \\
    """
}




workflow variant_call {
  take:
    sample_bams
    intervals
    panel_of_normals
    germline_resource
    refgenome
    args

  main:

    // best practices from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
    // How to Call somatic mutations using GATK4 Mutect2

    sample_bams_w_key = sample_bams.map{ it -> ["${it.getSimpleName()[0..-4]}", it]}
    bam_pairs = sample_bams_w_key.groupTuple(size: 2, sort: true).map{it -> it[1]}
    bam_pairs.view()

    refgenome_index = index_fasta(args.refgenome).fasta_index
    refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict

    indices = Channel.of(germline_resource, panel_of_normals)
    indices_mid = gatk_indexfeaturefile(indices).known_sites_index
    indices.view()
    indices_mid.map{it -> it.getSimpleName() }.view()
    println "${file(panel_of_normals).getSimpleName()}"


    // panel_of_normals_index = gatk_indexfeaturefile(panel_of_normals).known_sites_index
    panel_of_normals_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(panel_of_normals).getSimpleName()}" }
    germline_resource_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(germline_resource).getSimpleName()}" }

    mut = gatk_mutect(bam_pairs, intervals, 
	    panel_of_normals, panel_of_normals_index,
	    germline_resource,  germline_resource_index,  
	    refgenome, refgenome_index, refgenome_dict)
    om = gatk_learn_readorientationmodel(mut.f1r2_data).orientation_model

    // for tumors and for normals
    // todo, check that ew can use the same germline resource for pileup here
    all_pileups = gatk_getpileupsummaries(sample_bams, mut.vcf, intervals, germline_resource, germline_resource_index)

    //matched_pileups = all_pileups.groupTuple(size: 2, sort: true)

    //contamination = gatk_calculate_contamination(matched_pileups, segments_table)
    //filtered_vcf = gatk_filter_calls(mut.vcf, segments_table, contamination, om)
    filtered_vcf = []

  emit:
    vcf = filtered_vcf
}

workflow create_pon {
  take:
    args
  main:
    vcf = filtered_vcf
}
