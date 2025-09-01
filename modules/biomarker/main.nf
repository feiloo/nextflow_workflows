include { msisensor_pro } from "$NEXTFLOW_MODULES/biomarker/msisensor_pro.nf"
include { sequenza } from "$NEXTFLOW_MODULES/biomarker/sequenza_scarhrd.nf"
include { collect_hs_metrics } from "$NEXTFLOW_MODULES/biomarker/picard.nf"

workflow analyse_biomarkers {
  take:
    bam_pairs_w_idx
    refgenome
    refgenome_index
    refgenome_dict
    intervals
  main:
    out = msisensor_pro(bam_pairs_w_idx, refgenome)

    //seqq = sequenza(matched_bams, refgenome)

    bams = bam_pairs_w_idx.flatMap{normal_bam, normal_bai, tumor_bam, tumor_bai -> [normal_bam, tumor_bam]}
    //hs_metrics = collect_hs_metrics(bams,refgenome,refgenome_index,refgenome_dict, intervals).hs_metrics
    hs_metrics = Channel.empty()

    emit:
      matched_preproc_bams = out.matched_preproc_bams
      msi_csv = out.msi_csv
      hs_metrics = hs_metrics
      //indices = out.indices
}
