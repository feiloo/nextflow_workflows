nextflow.enable.dsl=2

def samplename_from_filename = { filename ->
    filename.split('-')[0]//.join('-')
}

def samplename_from_path = { path ->
    samplename_from_filename(path.name)
}

def year_from_path = { path ->
	file(path).name.split('-')[1]
}

process clc_workflow_batch {
  // imports the files into clc from an import/exportdir
  secret 'CLC_HOST'
  secret 'CLC_USER'
  secret 'CLC_PSW'

  container 'clc_client:latest'
  input:
    val files 
    val(clc_import_dir)
    val(clc_export_dir)
    val(clc_destdir)
    val(workflow_name)

  output:
    stdout emit: output
    val vcf_output, emit: vcf_output

  script:
  def arg_pref = "--workflow-input--5--select-files"
  def args2 = []
  files.each{f -> 
        args2 << "${arg_pref} \"clc://serverfile/${clc_import_dir}/${f[1].name}\"" 
	}

  def argstring = args2.join(" ")
  def destdir = clc_destdir
  def workflow_name = workflow_name

  def outdir = clc_export_dir
  def vcf_output = []

  files.each{f -> 
  	vcf_output << outdir + f[1].getSimpleName()[0..-8] + ' (paired) Unfiltered Variants-2.vcf.gz'
	}

  """
  #echo bla
  # clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${workflow_name} --workflow-input--5--import-command ngs_import_illumina -d ${destdir} ${argstring}
  # echo "-S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${workflow_name} --workflow-input--5--import-command ngs_import_illumina -d ${destdir} ${argstring}""
  """

}

process clc_workflow_single {
  // imports the files into clc from an import/exportdir
  secret 'CLC_HOST'
  secret 'CLC_USER'
  secret 'CLC_PSW'

  // use local container only for now
  container 'docker://localhost/clc_client:latest'

  input:
    tuple val(dna_read1), val(dna_read2), val(rna_read1), val(rna_read2) 
    val(clc_import_dir)
    val(clc_export_dir)
    val(clc_destdir)
    val(workflow_name)

  output:
    stdout emit: output

    script:

    def destdir = clc_destdir
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    def n3 = file(rna_read1).name
    def n4 = file(rna_read2).name
    def samplename = samplename_from_filename(n1)

    """
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A mkdir -t "${destdir}" -n "${samplename}"
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${workflow_name} \\
        --dna-reads-import-command ngs_import_illumina \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n1}\"  \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n2}\"  \\
        --rna-reads-import-algo-id ngs_import_illumina \\
        --rna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n3}\"  \\
        --rna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n4}\"  \\
        -d "${destdir}/${samplename}"
    """

    stub:

    def destdir = clc_destdir
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    def n3 = file(rna_read1).name
    def n4 = file(rna_read2).name
    def samplename = samplename_from_filename(n1)

    """
    echo "${destdir}"
    echo "${samplename}"
    echo "${n1}"
    echo "${n2}"
    echo "${n3}"
    echo "${n4}"
    echo "${destdir}/${samplename}"
    """
}

process copyfiles {
  // copies files into an clc import/exportdir
  input:
    val(read)
    val(nas_import_dir)
  output:
    val(ouf)

  script:
  def prefix = nas_import_dir
  def filen = file(read).getName()
  ouf = file(prefix)/filen

  """
  cp -n "${read}" "${ouf}"
  """

  stub:
  def prefix = nas_import_dir
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
  cp -r "${vcfin}" "${filen}"
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

def row_to_dict = { row ->
    [
    dna_read1: row[0],
    dna_read2: row[1],
    rna_read1: row[2],
    rna_read2: row[3]
    ]
}


workflow clc_nextflow {
  take:
    samplesheet
    clc_import_dir
    clc_export_dir
    clc_destdir
    clc_workflow_name
    nas_import_dir

  main:
    samples = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: true)

    // check that samplenames match and are in the
    def check_row = { row -> 
        def sid = file(row.dna_read1).name.split('-')[0]
	def year = year_from_path(file(row.dna_read1).name)
	def expected_row = [
		"${sid}-${year}-DNA_1.fq.gz",
		"${sid}-${year}-DNA_2.fq.gz",
		"${sid}-${year}-RNA_1.fq.gz",
		"${sid}-${year}-RNA_2.fq.gz",
	]
	def actual_row = [
		file(row.dna_read1).name,
		file(row.dna_read2).name,
		file(row.rna_read1).name,
		file(row.rna_read2).name,
	]
	if (actual_row != expected_row){
		throw new Exception("row doesnt match naming scheme ${row}")
	}
    }

    samples.subscribe{ row -> check_row(row) }

    sample_names = samples.map{ it -> [
    	"${it.dna_read1}", "${it.dna_read2}", 
	    "${it.rna_read1}", "${it.rna_read2}"]
	}

    sample_files = samples.flatMap{ it -> [
    	"${it.dna_read1}", "${it.dna_read2}", 
	    "${it.rna_read1}", "${it.rna_read2}"]
	}


    // use the files tuple for synchronizing the staging
    files = copyfiles(sample_files, nas_import_dir)
    staged_reads = files.map{it -> ["${samplename_from_path(it)}", it]}.groupTuple(by: 0, size: 4, sort: true).map{it -> it[1]}

    staged_reads.subscribe{ row -> check_row(row_to_dict(row)) }

    out = clc_workflow_single(staged_reads, 
    	clc_import_dir, clc_export_dir,
    	clc_destdir, clc_workflow_name)

    //reads = samplechannels.reads1.mix(samplechannels.reads2)
    //samples = files.groupTuple(size:2).buffer(size: 2)
    /*
    out = clc_workflow_batch(files.buffer(size: 2), args)
    vcfs = copyback(out.vcf_output.flatten().unique())
    sheet = writesamplesheet(vcfs.collect())
    sheet.view()

    emit:
      vcfs
      sheet
    */
}


workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  clc_nextflow(args.samplesheet, 
  	args.clc_import_dir, 
	args.clc_export_dir,
	args.clc_destdir,
	args.workflow_name,
	args.nas_import_dir
	)
}
