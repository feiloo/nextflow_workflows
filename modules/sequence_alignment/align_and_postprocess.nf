include { bwa_index_refgenome; bwa_align } from "$NEXTFLOW_MODULES/sequence_alignment/bwa.nf"
include { sam_to_bam; sort_bam; index_bam } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { publish } from "$NEXTFLOW_MODULES/sequence_alignment/utils.nf"

workflow bwa_tumor_only {
  take:
    args
  main:
    samplesheet = args.samplesheet
    // we require both samples to run the analysis
    // therefore 1 row requires/contains both, so its easier to read the samplesheet
    header = ['sample_id', 'tumor_read1', 'tumor_read2']

    csv_channel = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)
    all_reads = csv_channel.map{ it -> [it[0], it[1]] }

    bwa_idx = bwa_index_refgenome(args.refgenome)
    sams = bwa_align(all_reads, args.refgenome, bwa_idx.amb, bwa_idx.ann, bwa_idx.bwt, bwa_idx.pac, bwa_idx.sa)
    bams = sam_to_bam(sams)
    sorted_bams = sort_bam(bams)
    bam_indices = index_bam(sorted_bams)
    output_files = bam_indices.mix(sorted_bams)
    //publish(output_files, args)
}


workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  bwa_tumor_only(args)
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
