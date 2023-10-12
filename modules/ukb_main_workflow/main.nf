nextflow.enable.dsl=2

include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"
include { clc_nextflow } from "$NEXTFLOW_MODULES/clc_nextflow"
include { CIOABCD_VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation"
//include { SAREK } from "$NEXTFLOW_MODULES/sarek_wrapper"


workflow {
  //def args = params
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  args['input'] = '/opt/cio/samplesheet_testdata.csv'
  args['fasta'] = '/data/reference/ncbi_grch38_wholegenomefasta_genome.fa'
  args['vep_cache'] = '/data/reference/vep_cache/108/indexed/'


  //def workflow_name = args.workflow_name
  if(args.workflow_variation == 'arriba'){
  	arriba_nextflow(args)
  }
  else if(args.workflow_variation == 'clc'){
  	clc_nextflow(args)
  }
  else if(args.workflow_variation == 'sarek'){
  	SAREK(args)
  }
  else if(args.workflow_variation == 'variantinterpretation'){
  	CIOABCD_VARIANTINTERPRETATION(args)
  }
  else {
    println "invalid value for parameter workflow_variation"
    System.exit(1)
  }
}
