nextflow.enable.dsl=2

process arriba {
  container 'docker.io/uhrigs/arriba:2.4.0'
  cpus 8
  memory '50 GB'
  stageInMode 'copy'

  input: 
    tuple val(sampleid), path(read1), path(read2)
    path reference_data

  output:
    tuple val(sampleid), path('fusions.tsv'), path('fusions.discarded.tsv'), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')

  script:
  """
  touch sample_${sampleid}_info
  echo ${sampleid} ${read1} ${read2}
  /arriba_v2.4.0/run_arriba.sh \
  	${reference_data}/STAR_index_GRCh38_GENCODE38/ \
	${reference_data}/GENCODE38.gtf \
	${reference_data}/GRCh38.fa \
	${reference_data}/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \
	${reference_data}/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
  	${reference_data}/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
	${task.cpus} \
  	${read1} \
	${read2}
  rm -rf ${reference_data}
  """
}

process draw_fusions {
  container 'docker.io/uhrigs/arriba:2.4.0'
  cpus 1
  memory '4 GB'
  //stageInMode 'copy'

  input: 
    tuple val(sampleid), path('fusions.tsv'), path('fusions.discarded.tsv'), path('Aligned.sortedByCoord.out.bam'), path('Aligned.sortedByCoord.out.bam.bai')
    path reference_data

  output:
    tuple val(sampleid), path('fusions_report.pdf')

  script:
  """
  touch sample_${sampleid}_info
  /arriba_v2.4.0/draw_fusions.R \
  	--fusions=fusions.tsv \
	--annotation=${reference_data}/GENCODE38.gtf \
	--alignments=Aligned.sortedByCoord.out.bam \
	--cytobands=${reference_data}/cytobands_hg38_GRCh38_v2.4.0.tsv \
	--proteinDomains=${reference_data}/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
	--output=fusions_report.pdf
  """
}

process publish {
  cpus 1
  memory '1 GB'

  input:
    tuple val(sampleid), path('fusions.tsv'), path('fusions.discarded.tsv'), path('fusions_report.pdf')
    val(args)


  publishDir "${args.output_dir}", mode: 'copy', overwrite: false

  output:
    path 'sample_*'

  script:
  """
  mkdir sample_${sampleid}
  cp fusions.tsv fusions.discarded.tsv fusions_report.pdf sample_${sampleid}
  cd sample_${sampleid}
  echo ${sampleid} >> sample_id_info.txt
  """
    
}

workflow arriba_nextflow {
  take:
    args

  main:
    def samplesheet = args.samplesheet
    samples = Channel.fromPath(samplesheet).splitCsv(header: true)

    fusion_results = arriba(samples,args.reference_data)
    fusions_report = draw_fusions(fusion_results, args.reference_data)

    fusions_tables = fusion_results.map{t -> [t[0], t[1], t[2]]}
    results = fusions_tables.join(fusions_report)
    publish(results, args)


  //emit:
  //  fusion_results
}

workflow {
  def args = params
  arriba_nextflow(args)
}
