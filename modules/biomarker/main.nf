include { msisensor_pro } from "$NEXTFLOW_MODULES/biomarker/msisensor_pro.nf"
include { sequenza } from "$NEXTFLOW_MODULES/biomarker/sequenza_scarhrd.nf"
include { collect_hs_metrics } from "$NEXTFLOW_MODULES/biomarker/picard.nf"

workflow analyse_biomarkers {
  take:
    bams
    refgenome
    refgenome_index
    refgenome_dict
    intervals
    scar_hrd_header
    
  main:
    out = msisensor_pro(bams, refgenome)
    matched_bams = out.matched_preproc_bams

    seqz = sequenza(matched_bams, refgenome, scar_hrd_header)

    out2 = collect_hs_metrics(bams.flatten(),refgenome,refgenome_index,refgenome_dict, intervals)

    emit:
      matched_preproc_bams = out.matched_preproc_bams
      msi_csv = out.msi_csv
      hs_metrics = out2.hs_metrics
      indices = out.indices
      sequenza_results = seqz.sequenza_results
      hrd_results = seqz.scar_hrd_results
}
