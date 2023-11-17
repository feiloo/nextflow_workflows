nextflow.enable.dsl=2

process scan {
    input:
    path(refgenome)

    output:
    path("${refgenome}_msi.list"), emit:refgenome_microsatellites

    script:
    """
    msisensor-pro scan \\
    	-d ${refgenome} \\
	-o ${refgenome}_msi.list
    """
}

process eval_msi {
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
    args
  main:
    matched_sams
    sites = scan(refgenome).refgenome_microsatellites
    eval_msi(matched_sams, sites)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}
