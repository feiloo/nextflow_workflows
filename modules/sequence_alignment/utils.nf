process publish {
  memory '1 GB'

  //publishDir "${output_dir}/outputs", mode: 'copy', overwrite: false
  maxForks 1
  cache false

  input:
    path(inputfile)
    val(output_dir)
  output:
    path(inputfile)


  script:
  """
  mkdir -p ${output_dir}/outputs/
  # ensure that files are fully staged out, by rereading, and make it neary-atomic by writing to a tmp file and renaming

  echo publishing "${inputfile}" to "${output_dir}/outputs/${inputfile.getName()}"

  if [ -f "${output_dir}/outputs/${inputfile.getName()}" ]; then
      echo error, file already exists
      exit 1
  fi

  cp --no-clobber "${inputfile}" "${output_dir}/outputs/tmp.${inputfile.getName()}"
  sync
  # cmp is better than checksums because it fails earlier
  cmp ${inputfile} "${output_dir}/outputs/tmp.${inputfile.getName()}"
  sync
  mv --no-clobber "${output_dir}/outputs/tmp.${inputfile.getName()}" "${output_dir}/outputs/${inputfile.getName()}"
  sync
  echo successfully published "${inputfile}" to "${output_dir}/outputs/${inputfile.getName()}"
  """
}

process create_report {
  memory '1 GB'

  //publishDir "${output_dir}/outputs", mode: 'copy', overwrite: false
  maxForks 4
  cache false

  input:
    path("*.json")
  output:
    path("combined_report")


  script:
  """
  cat *.json >> combined_report
  """
}
