process publish {
  memory '1 GB'

  publishDir "${output_dir}/outputs", mode: 'copy', overwrite: false

  input:
    path(inputfile)
    val(output_dir)
  output:
    path(inputfile)


  script:
  """
  mkdir -p ${output_dir}/outputs/
  echo published ${inputfile} to ${output_dir}/outputs/${inputfile}
  """
}

