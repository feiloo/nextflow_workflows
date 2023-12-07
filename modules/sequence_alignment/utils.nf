process publish {
  // a process for explicit publishing of files
  // usefull for possible postprocessing of them later
  memory '1 GB'

  input:
    path(resultfile)
    val(args)
  output:
    path("out/${resultfile}")
    val(args)

  publishDir "${args.output_dir}", mode: 'copy', overwrite: false, enabled: true

  script:
  log.info "output file ${resultfile} published to output dir: ${args.output_dir}"
  """
  mkdir out
  ln -s ../"${resultfile}" "out/${resultfile}"
  """
}

