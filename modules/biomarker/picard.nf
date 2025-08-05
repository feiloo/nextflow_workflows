process collect_hs_metrics {
    conda "bioconda::picard=3.2.0"
    container "quay.io/biocontainers/picard:3.2.0"

    //memory "258 GB"
    // roughly 0.6 times filesize as ram consumption
    memory "${((task.attempt * 0.2 + 0.6) * bamfile.size() / 1000000000)} GB"

    input:
        path(bamfile)
        path(ref_gatk38_fasta)
        path(ref_gatk38_fai)
        path(ref_gatk38_dict)
        path(mappable_38_IntervalList)

    output:
        path("${bamfile.getSimpleName()}_bam_mertics.csv"), emit: hs_metrics

    script:
    """
    mkdir tmp
    picard BedToIntervalList -I ${mappable_38_IntervalList} \
	    -O list.interval_list -SD Homo_sapiens_assembly38.dict

    picard -Xmx256g CollectHsMetrics \
     -I ${bamfile} \
     -O ${bamfile.getSimpleName()}_bam_mertics.csv \
     -R ${ref_gatk38_fasta} \
     -TARGET_INTERVALS list.interval_list \
     -BAIT_INTERVALS list.interval_list \
     --TMP_DIR tmp
    """
}


workflow picard {
    
    take:
        bams
        ref_gatk38_fasta_ch
        ref_gatk38_fai_ch
        ref_gatk38_fai_dict
        mappable_38_IntervalList_ch

    main:
        out = collect_hs_metrics(bams,refgenome,refgenome_index, refgenome_dict, intervals)

    emit:
        hs_metrics = out.hs_metrics
        
}
    
/*
workflow {
	csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sample_id), file(row.bamfile) }
	ref_gatk38_fasta_ch = Channel.value(params.ref_gatk38_fasta)
	ref_gatk38_fai_ch = Channel.value(params.ref_gatk38_fai)
	mappable_38_IntervalList_ch = Channel.value(params.mappable_38_IntervalList)
    
    picard(csv_ch,ref_gatk38_fasta_ch,ref_gatk38_fai_ch,mappable_38_IntervalList_ch)    
}
*/
