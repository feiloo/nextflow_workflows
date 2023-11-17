nextflow.enable.dsl=2

include { index_bam; sort_bam } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"

process scan {
    container 'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0'
    conda 'bioconda::msisensor-pro=1.2.0'

    input:
    path(refgenome)

    output:
    path("${refgenome.getSimpleName()}_msi.list"), emit:refgenome_microsatellites

    script:
    """
    msisensor-pro scan \\
    	-d ${refgenome} \\
	-o ${refgenome.getSimpleName()}_msi.list
    """
}

process eval_msi {
    container 'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0'
    conda 'bioconda::msisensor-pro=1.2.0'

    // sorted and indexed bams
    input:
    tuple path(tumor_bam), path(normal_bam)
    path(refgenome_microsatellites)

    output:
    path("${normal_bam.getSimpleName()}_msi.csv")
    

    script:
    """
    msisensor-pro msi \\
        -d ${refgenome_microsatellites} \\
	-n ${normal_bam} \\
	-t ${tumor_bam} \\
	-o ${normal_bam.getSimpleName()}_msi.csv
    """
}


workflow msi_annotate {
  take:
    matched_bams
    refgenome
  main:
    sites = scan(refgenome).refgenome_microsatellites

    matched_bams.flatten().view()
    all_bams = matched_bams.flatten()
    sorted = sort_bam(all_bams)
    // add filename based key, group and remove key
    sorted_w_key = sorted.map{it -> [it.getSimpleName().split('_')[0], it]}
    sorted_w_key.view()
    //preproc_bams = index_bam(sorted).groupTuple().map{it -> it[1]}
    //eval_msi(preproc_bams, sites)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}
