include { msisensor_pro } from "$NEXTFLOW_MODULES/biomarker/msisensor_pro.nf"
include { sequenza } from "$NEXTFLOW_MODULES/biomarker/sequenza_scarhrd.nf"
include { collect_hs_metrics } from "$NEXTFLOW_MODULES/biomarker/picard.nf"

workflow analyse_biomarkers {
  take:
    bams
    bam_pairings
    refgenome
    refgenome_index
    refgenome_dict
    intervals
  main:
    out = msisensor_pro(bams, bam_pairings, refgenome)
    matched_bams = out.matched_preproc_bams

    //seqq = sequenza(matched_bams, refgenome)

    out2 = collect_hs_metrics(bams.flatten(),refgenome,refgenome_index,refgenome_dict, intervals)

    emit:
      matched_preproc_bams = out.matched_preproc_bams
      msi_csv = out.msi_csv
      hs_metrics = out2.hs_metrics
      indices = out.indices
}
