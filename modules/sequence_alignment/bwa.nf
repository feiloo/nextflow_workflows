process bwa_index_refgenome {
    conda "bioconda::bwa=0.7.17"
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    storeDir '/data2/bwa_indices'

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    //memory '120 GB'

    input:
    path(refgenome)

    output:

    // file endings are optionally prefix with an f so they start with an non-numerals
    // as a convention
    path("${refgenome}.amb"), emit: amb
    path("${refgenome}.ann"), emit: ann
    path("${refgenome}.bwt"), emit: bwt
    path("${refgenome}.pac"), emit: pac
    path("${refgenome}.sa"), emit: sa

    script:
    """
    bwa index ${refgenome}
    """
}

process bwa_align {
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    memory '56 GB'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path(refgenome)
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt")
    path("${refgenome}.pac")
    path("${refgenome}.sa")


    output:
    path("${read1.getSimpleName()}.sam")

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    #bwa mem -t $task.cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam
    bwa mem -t $n_cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam
    """
}

process sam_to_bam {
    conda "bioconda::samtools=1.17"
    container 'biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path(samfile)

    output:
    path("${samfile.getBaseName()}.bam")

    script:
    """
    samtools view ${samfile} --bam --threads 12 -o ${samfile.getBaseName()}.bam
    """
}

process publish {
  cpus 1
  memory '1 GB'
  

  input:
    path(bamfile)
    val(args)

  publishDir "${args.output_dir}", mode: 'copy', overwrite: false

  script:
  println "finished workflow, outputs are in: ${args.output_dir}"
  """
  touch ${bamfile}
  """

}

workflow bwa_tumor_only {
  take:
    args
  main:
    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)

    bwa_idx = bwa_index_refgenome(args.refgenome)
    sams = bwa_align(csv_channel, args.refgenome, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt, bwa_idx.pac, bwa_idx.sa)
    bams = sam_to_bam(sams)
    publish(sams, args)
}

/*

    i1 = Channel.fromPath("${args.refgenome}.amb")
    i2 = Channel.fromPath("${args.refgenome}.ann")
    i3 = Channel.fromPath("${args.refgenome}.bwt")
    i4 = Channel.fromPath("${args.refgenome}.pac")
    i5 = Channel.fromPath("${args.refgenome}.sa")

    bwa_index_refgenome(args.refgenome)
    bwa_align(csv_channel, args.refgenome, i1, i2, i3, i4, i5)

*/

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  bwa_tumor_only(args)

}

