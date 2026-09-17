include { CHECKBEDFILE		                        } from 'bedfile/checkbedfile/main'
include { TAGROI                                    } from 'subworkflows/vcf_roi_tagging/main'
include { BCFTOOLS_INDEX                            } from 'nf-core/bcftools/index/main'
include { SAMTOOLS_DICT                             } from 'nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX                            } from 'nf-core/samtools/faidx/main'
include { VCFPROC                                   } from 'subworkflows/process_vcf/main'
include { ENSEMBLVEP_FILTERVEP as TRANSCRIPT_FILTER } from 'nf-core/ensemblvep/filtervep/main'
include { ENSEMBLVEP_VEP                            } from 'nf-core/ensemblvep/vep/main'
include { TSV_CONVERSION                            } from 'subworkflows/tsv_conversion/main'
include { TMB_CALCULATE		                        } from 'tmbcalculation/main'
include { UKB_TOOL                                  } from 'UKB_tool/main'
include { UKB_TOOL_ONCOKB                           } from 'UKB_tool/main'

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

    // Check parameter combinations and give warnings
    if (!params.vep) log.warn("WARNING: You deactivated VEP-based annotation. Downstream processes are working properly only with VEP-annotated VCF file as input!")
    if (!params.tsv && params.report) error("ERROR: Needs to create TSV file for generating HTML report.")
    if (!params.tsv && params.calculate_tmb) error("ERROR: Need to create TSV file for calculating TMB.")
    if (!params.bedfile && params.tag_roi) error("ERROR: Need to specify bedfile for region-of-interest tagging.")
    if (!params.bedfile && params.calculate_tmb) error("ERROR: Need to specify bedfile for calculating TMB.")
    if (!params.read_depth && params.calculate_tmb) error("ERROR: Need to specify the read_depth FORMAT field for calculating TMB.")

    // Channels for UKB filter
    refseq_list                = params.refseq_list        ? Channel.value(params.refseq_list)                           : []
    variantDBi                 = params.variantDBi         ? Channel.value(params.variantDBi)                            : []

    // Index vcf and reference files

    // create tbi index for vcf
    BCFTOOLS_INDEX ( ch_samplesheet )
    vcf_tbi = ch_samplesheet.join(BCFTOOLS_INDEX.out.tbi)

    // create sequence dictionary and faidx index of reference FASTA
    fasta_ref = ch_fasta.map { ch_fasta -> ['ref', ch_fasta] }
    SAMTOOLS_DICT( fasta_ref )
    SAMTOOLS_FAIDX( fasta_ref, [[], []] )

    CHECKBEDFILE ( ch_bedfile )

    // ROI-tagging of VCF entries
    TAGROI (    ch_bedfile,
                vcf_tbi)
    tagroi_vcf=TAGROI.out.vcf_tbi

    // VCF filtering and normalization
    VCFPROC (
            tagroi_vcf,
            ch_fasta
    )


    // TODO: deleted vcf merging, might readd later but then cleaner
    proc_vcf=VCFPROC.out.vcf
        .map { meta, vcf -> tuple( meta, vcf, []) }

    ENSEMBLVEP_VEP( proc_vcf,
                    ch_vep_genome,
                    ch_vep_species,
                    ch_vep_cache_version,
                    ch_vep_cache,
                    fasta_ref,
                    ch_vep_extra_files)
    ch_vcf = ENSEMBLVEP_VEP.out.vcf

    // Filtering for transcripts
    if ( params.transcriptfilter || (params.transcriptlist!=[]) ) {
        TRANSCRIPT_FILTER(  ch_vcf,
                            ch_transcriptlist
        )
        ch_vcf_tf = TRANSCRIPT_FILTER.out.output
    } else {
        ch_vcf_tf = ch_vcf
    }

    // Use custom filters to tag variants and create subsets
    if (params.custom_filters) {
        PRESETS_FILTER_REPORT ( ch_vcf_tf,
                                ch_custom_filters)
        ch_vcf_tag = PRESETS_FILTER_REPORT.out.vcf
    } else {
        ch_vcf_tag = ch_vcf_tf
    }


    // MODULE: TSV conversion with vembrane table

    TSV_CONVERSION (ch_vcf_tag,
                    ch_annotation_fields
    )
    ch_tsv = TSV_CONVERSION.out.tsv

    // MODULE: TMB calculation
    /*
    somatic_files = TSV_CONVERSION.out.tsv.filter { meta, file ->
                         !meta.id.contains('_Tpavegermline')
                    }

    TMB_CALCULATE ( somatic_files,
                    ch_bedfile
                )
    */

    // MODULE: UKB filter and annotation

    ukb_results = Channel.empty()
    use_oncokb_token = use_proprietary

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

    emit:
    ukb_results

}
