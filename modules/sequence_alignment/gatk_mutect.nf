// implement best-practice somatic variant calling
// reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2

include { index_bam; index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

process gatk_mutect {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory '56 GB'

    input:
        tuple path(normal_bam), path(normal_bam_index), path(tumor_bam), path(tumor_bam_index)
	path(intervals)
	path(panel_of_normals)
	path(panel_of_normals_index)
	path(germline_resource)
	path(germline_resource_index)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      tuple path("${tumor_bam.getSimpleName()}_unfiltered.vcf"), path("${tumor_bam.getSimpleName()}_unfiltered.vcf.stats"), emit: vcf
      path("${tumor_bam.getSimpleName()}_f1r2_data.tar.gz"), emit: f1r2_data

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    def normal_samplename = "${normal_bam.getSimpleName()}"
    """
    mkdir -p tmp
    gatk Mutect2 \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
    	--native-pair-hmm-threads ${2*n_cpus} \\
	--tmp-dir tmp \\
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

    memory '56 GB'

    input:
        path(f1r2_data)

    output:
      path("${f1r2_data.getSimpleName()}_model.tar.gz"), emit: orientation_model

    script:
    """
    mkdir -p tmp
    gatk LearnReadOrientationModel \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--input "${f1r2_data}" \\
	--output "${f1r2_data.getSimpleName()}_model.tar.gz"
    """
}


process gatk_getpileupsummaries {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory '56 GB'

    input:
        tuple path(sample_bam), path(sample_bam_index)
        path(genomic_intervals)
	// prepared from like gnomAD with pop-freqs in the info field
        path(variant_frequency_vcf)
        path(variant_frequency_vcf_index)

    output:
      path("${sample_bam.getSimpleName()}_pileup.table"), emit: table

    script:
    """
    mkdir -p tmp
    gatk GetPileupSummaries \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	-I ${sample_bam} \\
	-L ${genomic_intervals} \\
	-V ${variant_frequency_vcf} \\
	-O "${sample_bam.getSimpleName()}_pileup.table"
    """
}

process gatk_calculate_contamination {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory '56 GB'

    input:
        tuple path(normal_pileups), path(tumor_pileups)

    output:
      tuple path("${tumor_pileups.getSimpleName()}_contamination.table"), path("${tumor_pileups.getSimpleName()}_segments.tsv"), emit: contamination_table_and_segments

    script:
    """
    mkdir -p tmp
    gatk CalculateContamination \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--input ${tumor_pileups} \\
	-matched ${normal_pileups} \\
	--tumor-segmentation ${tumor_pileups.getSimpleName()}_segments.tsv \\
	--output "${tumor_pileups.getSimpleName()}_contamination.table"
    """
}

process gatk_filter_calls {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    memory '56 GB'

    input:
        tuple path(sample_vcf), path(contamination_table), path(tumor_segments) path(orientation_model)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      path("${sample_vcf.getSimpleName()}_filtered.vcf"), emit: vcf

    script:
    """
    mkdir -p tmp
    gatk FilterMutectCalls \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--variant ${sample_vcf} \\
	--ob-priors ${orientation_model} \\
	--contamination-table-segmentation ${contamination_table} \\
	--tumor-segmentation ${tumor_segments} \\
	--output "${sample_vcf.getSimpleName()}_filtered.vcf" \\
	-R "${refgenome}"
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
    mkdir -p tmp
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

    /*
    bam_pairs_w_key = sample_bams_w_key.groupTuple(size: 2, sort: true)
    bam_pairs = bam_pairs_w_key.map{it -> it[1]}
    */
    //bam_pairs.view()

    sample_bams_w_key = sample_bams.map{ it -> ["${it.getSimpleName()}", it]}
    sample_bams_idx_w_key = index_bam(sample_bams).map{ it -> ["${it.getSimpleName()}", it] }
    sample_bams_w_indices = sample_bams_w_key.join(sample_bams_idx_w_key).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }
    bam_pairs = sample_bams_w_indices.groupTuple(by: 0, size:2, sort:{it[0]}).map{ it -> [it[1][0][0], it[1][0][1], it[1][1][0], it[1][1][1]] }
    bam_pairs.view()

    refgenome_index = index_fasta(args.refgenome).fasta_index
    refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict

    indices = Channel.of(germline_resource, panel_of_normals)
    indices_mid = gatk_indexfeaturefile(indices).known_sites_index

    // indices_mid.map{it -> it.getSimpleName() }.view()
    // println "${file(panel_of_normals).getSimpleName()}"


    // panel_of_normals_index = gatk_indexfeaturefile(panel_of_normals).known_sites_index
    panel_of_normals_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(panel_of_normals).getSimpleName()}" }
    germline_resource_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(germline_resource).getSimpleName()}" }


    mut = gatk_mutect(bam_pairs, intervals, 
	    panel_of_normals, panel_of_normals_index,
	    germline_resource,  germline_resource_index,  
	    refgenome, refgenome_index, refgenome_dict)

    vcf = mut.vcf.map{it -> it[0]}
    om = gatk_learn_readorientationmodel(mut.f1r2_data).orientation_model

    // for tumors and for normals
    // todo, check that ew can use the same germline resource for pileup here
    vcf.view()

    //sample_bams_w_small_key = sample_bams.map{it -> ["${it[0].getSimpleName().split('_')[0]}", it]}
    //bams_and_vcf = mut.vcf.map{it -> ["${it[0].getSimpleName().split('_')[0]}",it]}.combine(sample_bams_w_small_key, by: 0).map{it -> it[1..2]}
    //bams_and_vcf.view()

    sample_bams_w_indices_no_key = sample_bams_w_indices.map{it -> it[1]}
    all_pileups = gatk_getpileupsummaries(sample_bams_w_indices_no_key, intervals, germline_resource, germline_resource_index).table

    // group pileups by samplename
    matched_pileups = all_pileups.map{it -> ["${it.getSimpleName().split('_')[0]}", it]}.groupTuple(size: 2, sort: true).map{it -> it[1]}
    matched_pileups.view()

    contamination_table_and_segments = gatk_calculate_contamination(matched_pileups).contamination_table_and_segments
    c_w_key = contamination_table_and_segments.map{it -> ["{it.getSimpleName().split('_')[0]}", it]}

    vcf_w_key = vcf.map{it -> ["{it.getSimpleName().split('_')[0]}", it]}
    om_w_key = om.map{it -> ["{it.getSimpleName().split('_')[0]}", it]}
    vcf_w_filter_data = vcf_w_key.join(c_w_key).join(om_w_key).flatten().map{it -> it[1..4]}
    vcf_w_filter_data.view()
    

    filtered_vcf = gatk_filter_calls(vcf_w_filter_data, args.refgenome, refgenome_index, refgenome_dict).vcf

  emit:
    vcf = filtered_vcf
}

workflow create_pon {
  take:
    args
  main:
    vcf = filtered_vcf
}
