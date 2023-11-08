process publish {
  cpus 1
  memory '1 GB'
  

  input:
    path(resultfile)
    val(args)
  output:
    path(bamfile)

  publishDir "${args.output_dir}", mode: 'copy', overwrite: false

  script:
  println "finished workflow, outputs are in: ${args.output_dir}"
  """
  touch ${resultfile}
  """
}

