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
    memory "16 GB"

    // sorted and indexed bams
    input:
    tuple path(normal_bam), path(normal_bam_bai), path(tumor_bam), path(tumor_bam_bai)
    path(refgenome_microsatellites)

    output:
    path("${normal_bam.getSimpleName()}_${tumor_bam.getSimpleName()}_msi.csv"), emit: csv
    path("${normal_bam.getSimpleName()}_${tumor_bam.getSimpleName()}_msi.csv_dis")
    path("${normal_bam.getSimpleName()}_${tumor_bam.getSimpleName()}_msi.csv_germline")
    path("${normal_bam.getSimpleName()}_${tumor_bam.getSimpleName()}_msi.csv_somatic")
    

    script:
    """
    msisensor-pro msi \\
        -d ${refgenome_microsatellites} \\
	-n ${normal_bam} \\
	-t ${tumor_bam} \\
	-o ${normal_bam.getSimpleName()}_${tumor_bam.getSimpleName()}_msi.csv
    """
}


workflow msisensor_pro {
  take:
    bam_pairs_w_idx
    refgenome
  main:
    sites = scan(refgenome).refgenome_microsatellites

    /*
    all_bams = bams.flatten()
    sorted = sort_bam(all_bams)

    indices = index_bam(sorted)
    // add filename based key, group and remove key
    preproc_bams = sorted.mix(indices).map{it -> [it.getSimpleName().split('_')[0], it]}
    // group and sort by filename 
    matched_preproc_bams = preproc_bams.groupTuple(size: 4).map{it -> it[1].sort{ e1, e2 -> e1.getName() <=> e2.getName() }}

    // assert that the pattern order and names match:
    // _N_1.bam, _N_1.bam.bai, _T_1.bam, _T_1.bam.bai
    matched_preproc_bams.subscribe{ it -> 
    	assert it.join(",") ==~ /.*_N_.*_1\.bam.*_N_.*_1\.bam\.bai.*_T_.*_1\.bam.*_T_.*_1\.bam\.bai/
	}

    sample_bams = sorted

    // key is XXX_T_[FFPE|BLOOD|FF]_1
    sample_bams_w_key = sample_bams.map{ it -> ["${it.getSimpleName()}", it]}
    sample_bams_idx_w_key = index_bam(sample_bams).map{ it -> ["${it.getSimpleName()}", it] }
    // join bams and bai files
    sample_bams_w_indices_w_key = sample_bams_w_key.join(sample_bams_idx_w_key)
    sample_bams_w_indices_w_key.view()
    sample_bams_w_indices = sample_bams_w_indices_w_key.map{it -> it[1..-1]}
    
    indices = sample_bams_w_indices.flatMap{ bam, bai -> [bai] }.unique()


    bams_split = sample_bams_w_indices_w_key.branch { key, bam, bai ->
	// Match pattern: <alphanum>-<2-digit>_<N or T>_
	def normal = key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_N_.+/
	def tumor  = key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_T_.+/

	// Output channels
	normal: normal
	tumor:  tumor
    }

    // Join with normal bams: match normal_key == key
    with_normal = bam_pairings
	.join(bams_split.normal)  // emits [normal_key, tumor_key, normal_bam, normal_bai]
	.map { normal_key, tumor_key, normal_bam, normal_bai ->
	    [tumor_key, normal_bam, normal_bai]  // prepare for tumor join
	}

    // Join with tumor bams: match tumor_key == key
    bam_pairs_w_idx = with_normal
	.join(bams_split.tumor)  // finds matching tumor
	.map { tumor_key, normal_bam, normal_bai, tumor_bam, tumor_bai ->
	    [normal_bam, normal_bai, tumor_bam, tumor_bai]
	}
    */


    msi_csv = eval_msi(bam_pairs_w_idx, sites).csv

    matched_preproc_bams = bam_pairs_w_idx

    emit:
      matched_preproc_bams
      //indices
      msi_csv
}


workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  msisensor_pro(args.bams, args.refgenome)

}
