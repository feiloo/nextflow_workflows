nextflow.enable.dsl=2
/*
 * perform picard CollectHsMetrics for WGS Pilot
 */

process collect_hs_mertics {
    
    tag "${sample_id}"
    debug true

    publishDir "${outdir}/QC/${sample_id}", mode: "copy"

    input:
        tuple val(sample_id), path(bamfile)
        path(ref_gatk38_fasta)
        path(ref_gatk38_fai)
        path(mappable_38_IntervalList)

    output:
        tuple val(sample_id), path("${sample_id}_bam_mertics.csv"), emit: bamqc

    script:
    """
    java -jar -Xmx256g /home/pbasitta/miniconda3/envs/picard_wgs_pilot/share/picard-3.2.0-0/picard.jar CollectHsMetrics \
     -I ${bamfile} \
     -O ${sample_id}_bam_mertics.csv \
     -R ${ref_gatk38_fasta} \
     -TARGET_INTERVALS ${mappable_38_IntervalList} \
     -BAIT_INTERVALS ${mappable_38_IntervalList} \
     --TMP_DIR ${tmp_dir}  
    """
}

csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sample_id), file(row.bamfile) }
ref_gatk38_fasta_ch = Channel.value(params.ref_gatk38_fasta)
ref_gatk38_fai_ch = Channel.value(params.ref_gatk38_fai)
mappable_38_IntervalList_ch = Channel.value(params.mappable_38_IntervalList)
tmp_dir = params.tmp_dir
outdir = params.outdir

workflow picard {
    
    take:
        csv_ch
        ref_gatk38_fasta_ch
        ref_gatk38_fai_ch
        mappable_38_IntervalList_ch

    main:
        collect_hs_mertics(csv_ch,ref_gatk38_fasta_ch,ref_gatk38_fai_ch,mappable_38_IntervalList_ch)
        
}
    
workflow {
    
    picard(csv_ch,ref_gatk38_fasta_ch,ref_gatk38_fai_ch,mappable_38_IntervalList_ch)    
}
