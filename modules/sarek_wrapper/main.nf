include {NFCORE_SAREK} from './sarek' //addParams(input_restart: null)

workflow SAREK {
  take:
    args
  main:
    NFCORE_SAREK(args)
}

/*
workflow {
  main:
    def args = params
    SAREK(args)
}
*/
