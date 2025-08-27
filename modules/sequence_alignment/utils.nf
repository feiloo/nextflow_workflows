process publish {
  memory '1 GB'

  //publishDir "${output_dir}/outputs", mode: 'copy', overwrite: false
  maxForks 4
  cache false

  input:
    path(inputfile)
    val(output_dir)
  output:
    path(inputfile)


  script:
  """
  mkdir -p ${output_dir}/outputs/

  echo published ${inputfile} to ${output_dir}/outputs/${inputfile.getName()}
  cp --no-clobber ${inputfile} ${output_dir}/outputs/${inputfile.getName()}
  sha256sum ${inputfile} >> input_hash
  sha256sum ${output_dir}/outputs/${inputfile.getName()} >> output_hash
  """
}

