nextflow.enable.dsl=2


process clc_workflow_batch {
  secret 'CLC_HOST'
  secret 'CLC_USER'
  secret 'CLC_PSW'

  container 'clcclient:latest'
  input:
    val files 
    val args

  output:
    stdout emit: output
    val vcf_output, emit: vcf_output

  script:
  def arg_pref = "--workflow-input--5--select-files"
  def args2 = []
  files.each{f -> 
        args2 << "${arg_pref} \"clc://serverfile/${args.clc_import_dir}/${f[1].name}\"" 
	}

  def argstring = args2.join(" ")
  def destdir = args.clc_destdir
  def workflow_name = args.workflow_name

  def outdir = args.clc_export_dir
  def vcf_output = []

  files.each{f -> 
  	vcf_output << outdir + f[1].getSimpleName()[0..-8] + ' (paired) Unfiltered Variants-2.vcf.gz'
	}

  """
  echo bla
  # clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${workflow_name} --workflow-input--5--import-command ngs_import_illumina -d ${destdir} ${argstring}
  """

  stub:
  def prefix = args.nas_import_dir

  vcf_output = []
  files.each{f -> 
        vcf_output << "${prefix}/${f[1].name}.vcf"
	}

  """
  echo bla
  """
}


process copyfiles {
  input:
    tuple val(sid), val(read)
    val args
  output:
    tuple val(sid), val(ouf)

  script:
  def prefix = args.nas_import_dir
  def filen = file(read).getName()
  ouf = file(prefix)/filen

  """
  cp "${read}" "${ouf}"
  """

  stub:
  def prefix = args.nas_import_dir
  def filen = file(read).getName()
  ouf = file(prefix)/filen

  """
  touch "${ouf}"
  """
}

process copyback {
  input:
    val vcfin
  output:
    path filen

  script:

  filen = file(vcfin).name
  """
  cp "${vcfin}" "${filen}"
  """
  stub:
  filen = file(vcfin).name
  """
  touch "${filen}"
  """
}


process writesamplesheet {
  input:
    path vcfin
  output:
    path fil

  exec:

  fil = file(task.workDir + '/samplesheet.csv')
  fil.append("sample,vcf")

  vcfin.each{f ->
  	fil << f.name+','+f
	println f.name
  }
}


workflow clc_nextflow {
  take:
    args

  main:
    def samplesheet = args.samplesheet
    samples = Channel.fromPath(samplesheet).splitCsv(header: true)

    samples.multiMap{ row ->
          sampleid: "${row.sampleid}"
          reads1: ["${row.sampleid}", "${row.read1}"]
          reads2: ["${row.sampleid}", "${row.read2}"]
          }.set{samplechannels}

    reads = samplechannels.reads1.mix(samplechannels.reads2)
    reads.view()
    files = copyfiles(reads, args)
    //samples = files.groupTuple(size:2).buffer(size: 2)

    //samples.view()
    out = clc_workflow_batch(files.buffer(size: 2), args)
    vcfs = copyback(out.vcf_output.flatten().unique())
    sheet = writesamplesheet(vcfs.collect())
    sheet.view()

    emit:
      vcfs
      sheet
}


workflow {
  def args = params
  clc_nextflow(args)
}
