process bwamem2_index_refgenome {
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    // needs 28N GB, where N is the size of the uncompressed refseq in GB
    memory '120 GB'

    input:
    path(refgenome)

    output:

    // file endings are optionally prefix with an f so they start with an non-numerals
    // as a convention
    path("${refgenome}.0123"), emit: f0123
    path("${refgenome}.amb"), emit: amb
    path("${refgenome}.ann"), emit: ann
    path("${refgenome}.bwt.2bit.64"), emit: bwt_2bit_64
    path("${refgenome}.pac"), emit: pac

    script:
    """
    bwa-mem2 index ${refgenome}
    """
}

process bwamem2_align {
    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    memory '120 GB'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path(refgenome)
    path("${refgenome}.0123")
    path("${refgenome}.amb")
    path("${refgenome}.ann")
    path("${refgenome}.bwt.2bit.64")
    path("${refgenome}.pac")


    //output:
    path("${read1.getSimpleName()}.sam")

    script:
    """
    bwa-mem2 mem -t $task.cpus ${refgenome} ${read1} ${read2} -o ${read1.getSimpleName()}.sam
    """
}

workflow bwamem2_tumor_only {
  take:
    args
  main:
    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)

    i1 = Channel.fromPath("${args.refgenome}.0123")
    i2 = Channel.fromPath("${args.refgenome}.bwt.2bit.64")
    i3 = Channel.fromPath("${args.refgenome}.alt")
    i4 = Channel.fromPath("${args.refgenome}.ann")
    i5 = Channel.fromPath("${args.refgenome}.pac")

    bwamem2_align(csv_channel, args.refgenome, i1, i2, i3, i4, i5)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  bwamem2_tumor_only(args)

}

