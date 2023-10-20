nextflow.enable.dsl=2

include { fastqc } from "$NEXTFLOW_MODULES/sequence_alignment/fastqc.nf"

workflow sequence_alignment {
  take:
    args
  main:
    samplesheet = args.samplesheet
    header = ['sample_id', 'normal_read1', 'normal_read2', 'tumor_read1', 'tumor_read2']

    samples = Channel.fromPath(samplesheet).splitCsv(header: header, skip: 1)

    reads = samples.flatMap{
	row -> [ [row.sampleid, row.normal_read1],
		[row.sampleid, row.normal_read2],
		[row.sampleid, row.tumor_read1],
		[row.sampleid, row.tumor_read2]
		]
	}

    fastqc(reads)
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  sequence_alignment(args)

}
