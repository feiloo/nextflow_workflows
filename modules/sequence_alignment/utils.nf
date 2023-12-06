process publish {
  // a process for explicit publishing of files
  // usefull for possible postprocessing of them later
  memory '1 GB'

  input:
    path(resultfile)
    val(args)
  output:
    path("out/${resultfile}")

  publishDir "${args.output_dir}", mode: 'copy', overwrite: false

  script:
  println "finished workflow, outputs are in: ${args.output_dir}"
  """
  mkdir out
  ln -s "out/${resultsfile}" ${resultsfile}
  """
}

