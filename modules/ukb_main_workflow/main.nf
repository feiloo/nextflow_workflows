nextflow.enable.dsl=2

include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"
include { clc_nextflow } from "$NEXTFLOW_MODULES/clc_nextflow"
include { CIOABCD_VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation"
include { sequence_alignment } from "$NEXTFLOW_MODULES/sequence_alignment"
//include { SAREK } from "$NEXTFLOW_MODULES/sarek_wrapper"


workflow {
  //def args = params
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  println args

  if(args.workflow_variation == 'sequence_alignment'){
  	sequence_alignment(args)
  }
  else if(args.workflow_variation == 'arriba'){
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
    println "invalid value for parameter workflow_variation " + args.workflow_variation
    System.exit(1)
  }
}
