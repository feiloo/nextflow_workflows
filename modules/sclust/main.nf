nextflow.enable.dsl=2

process sclust_extract_read_ratio {
  cpus 1
  memory '4 GB'

  input: 
    tuple val(sampleid), path(bam_n), path(bam_n_index), path(bam_t), path(bam_t_index), val(chromosome)

  output:
    tuple val(sampleid), path("1_${chromosome}_bamprocess_data.txt"), emit: read_ratios

  script:
  """
  echo ${bam_n_index} ${bam_t_index}
  sclust bamprocess -n ${bam_n} -t ${bam_t} -part 1 -build hg19 -r "${chromosome}" -o 1
  """
}

process sclust_merge_chromosomes {
  cpus 1
  memory '4 GB'

  input: 
    tuple val(sampleid), path("1_*_bamprocess_data.txt")

  output:
    tuple val(sampleid), path("1_rcount.txt"), path("1_snps.txt"), emit: counts

  script:
  """
  sclust bamprocess -i 1 -o 1
  """
}


process sclust_copy_number_analysis {
  cpus 1
  memory '4 GB'

  input: 
    tuple val(sampleid), path('sample_mutations.vcf'), path("1_rcount.txt"), path("1_snps.txt") 

  output:
    tuple val(sampleid)

  script:
  """
  sclust cn -rc ${1_rcount.txt} -snp ${1_snps.txt} -vcf ${sample_mutations.vcf} -o 1
  """
}



workflow sclust_nextflow {
  take:
    args
    chromosomes

  main:
    def samplesheet = args.samplesheet
    // def header = ['sample_id', 'vcf', 'bam_n', 'bam_n_index', 'bam_t', 'bam_t_index']
    samples = Channel.fromPath(samplesheet).splitCsv(skip: 1)
    combs = samples.combine(chromosomes)

    results = sclust_extract_read_ratio(combs).read_ratios
    sample_chromosomes = results.groupTuple()
    sample_chromosomes.view()
    counts = sclust_merge_chromosomes(sample_chromosomes).counts

    sclust_copy_number_analysis(counts)



  //emit:
  //  fusion_results
}

workflow {
  def args = params
  def chromosomes = []
  for(int i=1;i<23;i++) {
    chromosomes << "chr$i"
  }
  chch = Channel.from(chromosomes)

  sclust_nextflow(args, chch)
}
