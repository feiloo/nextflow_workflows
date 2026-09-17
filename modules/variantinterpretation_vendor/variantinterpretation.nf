include { CHECKBEDFILE		                        } from '../modules/local/bedfile/checkbedfile/main'
include { TAGROI                                    } from '../subworkflows/local/vcf_roi_tagging/main'
include { BCFTOOLS_INDEX                            } from '../modules/nf-core/bcftools/index/main'
include { SAMTOOLS_DICT                             } from '../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX                            } from '../modules/nf-core/samtools/faidx/main'
include { CHECKVCF                                  } from '../subworkflows/local/check_vcf/main'
include { VCFPROC                                   } from '../subworkflows/local/process_vcf/main'
include { MERGE_VCFS                                } from '../subworkflows/local/merge_vcfs/main'
include { ENSEMBLVEP_FILTERVEP as TRANSCRIPT_FILTER } from '../modules/nf-core/ensemblvep/filtervep/main'
include { ENSEMBLVEP_VEP                            } from '../modules/nf-core/ensemblvep/vep/main'
include { TSV_CONVERSION                            } from '../subworkflows/local/tsv_conversion/main'
include { VARIANTFILTER as PRESETS_FILTER_REPORT    } from '../subworkflows/local/variantfilter/main'
include { HTML_REPORT                               } from '../subworkflows/local/html_report/main'
include { TMB_CALCULATE		                    } from '../modules/local/tmbcalculation/main'
include { UKB_FILTER                                } from '../modules/local/UKB_filter/main'
include { UKB_TOOL                                } from '../modules/local/UKB_tool/main'
include { UKB_TOOL_ONCOKB                                } from '../modules/local/UKB_tool/main'
include { ONCOKB_ANNOTATOR_UKB                      } from '../modules/local/oncokb_annotator_ukb/main'
include { WXS_ANNOTATION_UKB                        } from '../modules/local/wxs_annotation_ukb/main'


workflow VARIANTINTERPRETATION {

    take:
    ch_samplesheet
    ch_fasta
    ch_vep_cache
    ch_vep_cache_version
    ch_vep_genome
    ch_vep_species
    ch_vep_extra_files
    ch_annotation_fields
    ch_transcriptlist
    ch_datavzrd_config
    ch_annotation_colinfo
    ch_bedfile
    ch_custom_filters
    ch_library_type
    use_proprietary_arg

    main:
    // gather versions of each process
    ch_versions = Channel.empty()
    // gather QC reports for multiQC
    ch_multiqc_files = Channel.empty()
    // gather warnings
    ch_warnings = Channel.empty()

    // yikes, nextflow type coerces the argument values (true and false, ...) into boolean instead of keeping them strings
	if (use_proprietary_arg == null) {
      use_proprietary = false
	} else if (use_proprietary_arg == true) {
	  use_proprietary = true
	} else if (use_proprietary_arg == false) {
	  use_proprietary = false
	} else {
	  throw new IllegalArgumentException("use_proprietary must be boolean: true or false")
	}

	if (use_proprietary){
		  println 'enabeling proprietary features'
	}

    //
    // Check parameter combinations and give warnings
    //
    if (!params.vep) log.warn("WARNING: You deactivated VEP-based annotation. Downstream processes are working properly only with VEP-annotated VCF file as input!")
    if (!params.tsv && params.report) error("ERROR: Needs to create TSV file for generating HTML report.")
    if (!params.tsv && params.calculate_tmb) error("ERROR: Need to create TSV file for calculating TMB.")
    if (!params.bedfile && params.tag_roi) error("ERROR: Need to specify bedfile for region-of-interest tagging.")
    if (!params.bedfile && params.calculate_tmb) error("ERROR: Need to specify bedfile for calculating TMB.")
    if (!params.read_depth && params.calculate_tmb) error("ERROR: Need to specify the read_depth FORMAT field for calculating TMB.")

    // Channels for UKB filter
    refseq_list                = params.refseq_list        ? Channel.value(params.refseq_list)                           : []
    variantDBi                 = params.variantDBi         ? Channel.value(params.variantDBi)                            : []

    //
    // Index vcf and reference files
    //

    // create tbi index for vcf
    BCFTOOLS_INDEX ( ch_samplesheet )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
    vcf_tbi = ch_samplesheet.join(BCFTOOLS_INDEX.out.tbi)

    // create sequence dictionary and faidx index of reference FASTA
    fasta_ref = ch_fasta.map { ch_fasta -> ['ref', ch_fasta] }
    SAMTOOLS_DICT( fasta_ref )
    ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
    SAMTOOLS_FAIDX( fasta_ref, [[], []] )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    CHECKBEDFILE ( ch_bedfile )

    //
    // ROI-tagging of VCF entries
    //
    if (params.tag_roi && CHECKBEDFILE.out.bed_valid) {
        TAGROI (    ch_bedfile,
                    vcf_tbi)
        ch_versions = ch_versions.mix(TAGROI.out.versions)
        tagroi_vcf=TAGROI.out.vcf_tbi
    } else {
        tagroi_vcf=vcf_tbi
    }

    //
    // VCF filtering and normalization
    //
    VCFPROC (
            tagroi_vcf,
            ch_fasta
    )
    ch_versions = ch_versions.mix(VCFPROC.out.versions)

    //
    // Merging VCF files by groups
    //

    if (params.merge_vcfs) {
        MERGE_VCFS (
            VCFPROC.out.vcf
        )
        ch_versions = ch_versions.mix(VCFPROC.out.versions)

        proc_vcf=MERGE_VCFS.out.vcf
    } else {
        proc_vcf=VCFPROC.out.vcf
    }

    //
    // MODULE: VEP annotation
    //

    proc_vcf=proc_vcf
        .map { meta, vcf -> tuple( meta, vcf, []) }

    if (params.vep) {
        ENSEMBLVEP_VEP( proc_vcf,
                        ch_vep_genome,
                        ch_vep_species,
                        ch_vep_cache_version,
                        ch_vep_cache,
                        fasta_ref,
                        ch_vep_extra_files)
        ch_vcf = ENSEMBLVEP_VEP.out.vcf
        ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ENSEMBLVEP_VEP.out.report)
    } else {
        ch_vcf = proc_vcf
    }

    // Filtering for transcripts
    if ( params.transcriptfilter || (params.transcriptlist!=[]) ) {
        TRANSCRIPT_FILTER(  ch_vcf,
                            ch_transcriptlist
        )
        ch_vcf_tf = TRANSCRIPT_FILTER.out.output
        ch_versions = ch_versions.mix(TRANSCRIPT_FILTER.out.versions)
    } else {
        ch_vcf_tf = ch_vcf
    }

    // Use custom filters to tag variants and create subsets
    if (params.custom_filters) {
        PRESETS_FILTER_REPORT ( ch_vcf_tf,
                                ch_custom_filters)
        ch_vcf_tag = PRESETS_FILTER_REPORT.out.vcf
        ch_versions = ch_versions.mix(PRESETS_FILTER_REPORT.out.versions)
    } else {
        ch_vcf_tag = ch_vcf_tf
    }


    //
    // MODULE: TSV conversion with vembrane table
    //

    ukb_results = Channel.empty()

    TSV_CONVERSION (ch_vcf_tag,
                    ch_annotation_fields
    )
    ch_tsv = TSV_CONVERSION.out.tsv
    ch_versions = ch_versions.mix(TSV_CONVERSION.out.versions)

    //
    // MODULE: TMB calculation
    //
    somatic_files = TSV_CONVERSION.out.tsv.filter { meta, file ->
                         !meta.id.contains('_Tpavegermline')
                    }

    if ( params.bedfile && params.calculate_tmb ) {
            if ( CHECKBEDFILE.out.bed_valid ) {
                    TMB_CALCULATE ( somatic_files,
                                    ch_bedfile
                )
                ch_versions = ch_versions.mix(TMB_CALCULATE.out.versions)
        }
    }

    // MODULE: UKB filter

    def use_old_filter = false
    use_oncokb_token = use_proprietary

    if( use_old_filter == true ) {
        UKB_FILTER(ch_tsv, refseq_list, variantDBi, ch_library_type)
        ch_versions = ch_versions.mix(UKB_FILTER.out.versions)
        ch_filtered_variants = UKB_FILTER.out.variants_filtered_maf
        tmb = UKB_FILTER.out.tmb.map{it -> it[1]}
        ONCOKB_ANNOTATOR_UKB(ch_filtered_variants)
        annotated_variants = WXS_ANNOTATION_UKB(ONCOKB_ANNOTATOR_UKB.out.oncokb_out).annotated_variants
        ukb_results = ukb_results.mix(annotated_variants.map{ it -> it[1] } ).mix(tmb)
    } else {

        println "using new ukb_tool"
        //filtout = UKB_TOOL(ch_tsv, refseq_list, variantDBi, ch_library_type)
        if( use_oncokb_token == true ){
            println 'using onkokb token'
            filtout = UKB_TOOL_ONCOKB(ch_tsv, ch_library_type)
        }
        else {
            filtout = UKB_TOOL(ch_tsv, ch_library_type)
            }
        tmb = filtout.tmb.map{it -> it[1]}
        annotated_variants = filtout.annotated_variants
        ukb_results = ukb_results.mix(annotated_variants.map{ it -> it[1]} ).mix(tmb)
    }

    }

    emit:
    ch_versions
    ch_multiqc_files
    ch_warnings
    ukb_results

}
