include {NFCORE_SAREK} from './sarek'

workflow SAREK_WRAPPER {
  take:
    args
  main:
    println params.genome
    println params.genomes['GATK.GRCh38']
    NFCORE_SAREK(args)
}

workflow {
  main:
    def args = [:]
    for (param in params) { args[param.key] = param.value }
    for (param in params.genomes[params.genome]) { args[param.key] = param.value }

    SAREK_WRAPPER(args)
}
