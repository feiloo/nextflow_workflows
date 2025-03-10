nextflow.enable.dsl=2

include { index_bam; sort_bam } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"

process scan {
    container 'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0'
    conda 'bioconda::msisensor-pro=1.2.0'
    memory "150 GB"

    storeDir "$NEXTFLOW_STOREDIR/msisensorpro_scan"

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
    memory "10 GB"

    // sorted and indexed bams
    input:
    tuple path(normal_bam), path(normal_bam_bai), path(tumor_bam), path(tumor_bam_bai)
    path(refgenome_microsatellites)

    output:
    path("${normal_bam.getSimpleName()}_msi.csv"), emit: csv
    path("${normal_bam.getSimpleName()}_msi.csv_dis")
    path("${normal_bam.getSimpleName()}_msi.csv_germline")
    path("${normal_bam.getSimpleName()}_msi.csv_somatic")
    

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

    all_bams = matched_bams.flatten()
    sorted = sort_bam(all_bams)
    indices = index_bam(sorted)
    // add filename based key, group and remove key
    preproc_bams = sorted.mix(indices).map{it -> [it.getSimpleName().split('_')[0], it]}
    // group and sort by filename 
    matched_preproc_bams = preproc_bams.groupTuple(size: 4).map{it -> it[1].sort{ e1, e2 -> e1.getName() <=> e2.getName() }}

    // assert that the pattern order and names match:
    // _normal.bam, _normal.bam.bai, _tumor.bam, _tumor.bam.bai
    matched_preproc_bams.subscribe{ it -> 
    	assert it.join(",") ==~ /.*_normal\.bam.*_normal\.bam\.bai.*_tumor\.bam.*_tumor\.bam\.bai/
	}

    msi_csv = eval_msi(matched_preproc_bams, sites).csv
    emit:
      matched_preproc_bams
      msi_csv

}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}
