nextflow.enable.dsl=2


// example 10000-24-DNA_1.fq.gz -> 10000-24
// todo robustness/parsing
def samplename_from_filename = { filename ->
    filename.split('-')[0..1].join('-')
}

def samplename_from_path = { path ->
    samplename_from_filename(path.name)
}

def year_from_path = { path ->
    file(path).name.split('-')[1]
}


process clc_pancancer_dna_only {
  // imports the files into clc from an import/exportdir
  secret 'CLC_HOST'
  secret 'CLC_USER'
  secret 'CLC_PSW'

  time '80h'

  // use local container only for now
  container 'docker://localhost/clc_client:latest'

  input:
    tuple val(dna_read1), val(dna_read2)
    // windows path, where the windows-hosted clc server expects to import files from
    // must be the full mounted path if its on a network drive and properly slash escaped
    // \\host.domain\MOUNTPOINT\SHARENAME\IMPORT\...
    // is therefore prefixed wth clc://serverfile
    val(clc_import_dir)
    // windows path, where the windows-hosted clc server expects to export files to
    // same rules as for clc_import_dir
    // can be the same as clc_import_dir
    val(clc_export_dir)
    // clc path where clc should write workflow outputs (clc storage, e.g. clc://server/...)
    val(clc_destdir)
    // string that specifies the clc workflow to call as required by the clc command line tools
    val(clc_workflow_name)
    // unix path, where the network share is mounted, to copy the results from there
    val(nas_import_dir)
    val(nas_export_dir)
    // workflow run id
    val(workflow_run_id)

  output:
    stdout emit: output
    val(vcf_path), emit: vcf_path
    val(csv_path), emit: csv_path

    script:

    def run_id_clc = workflow_run_id.replace(":","_")
    def workflow_run_destdir = "${clc_destdir + '/' + run_id_clc}"
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    //todo add checks for filen naming scheme, n1, ... n4, should have the same prefix
    // this is used below as the a directory to group the outputs related to the sample
    def samplename = samplename_from_filename(n1)

    // compose name of clc output files (this depends on the clc workflow settings)
    def output_filename = "${samplename}"

    vcf_path = "${nas_export_dir}/${samplename}/${output_filename}.vcf"
    csv_path = "${nas_export_dir}/${samplename}/${output_filename}.csv"

    // clc requires to specify output folders, iirc doesnt allow setting filename with the cli tools
    """
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A mkdir -t "${clc_destdir}" -n "${workflow_run_id}"
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A mkdir -t "${workflow_run_destdir}" -n "${samplename}"
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${clc_workflow_name} \\
        --dna-reads-import-command ngs_import_illumina \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n1}\"  \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n2}\"  \\
        --export-table-csv-csv-destination \"clc://serverfile/${clc_export_dir}/${samplename}\"  \\
        --export-table-csv-custom-file-name "${samplename}.csv\"  \\
        --export-vcf-vcf-destination \"clc://serverfile/${clc_export_dir}/${samplename}\"  \\
        --export-vcf-custom-file-name "${samplename}.vcf\"  \\
        -d "${workflow_run_destdir}/${samplename}"
    """

    stub:

    def destdir = clc_destdir
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    def samplename = samplename_from_filename(n1)
    vcf_path = "${nas_import_dir}/${samplename}/${samplename}.vcf"
    csv_path = "${nas_import_dir}/${samplename}/${samplename}.csv"

    """
    echo "${destdir}"
    echo "${samplename}"
    echo "${n1}"
    echo "${n2}"
    echo "${destdir}/${samplename}"
    """
}


process clc_pancancer_dna_rna {
  // imports the files into clc from an import/exportdir
  secret 'CLC_HOST'
  secret 'CLC_USER'
  secret 'CLC_PSW'

  time '80h'

  // use local container only for now
  container 'docker://localhost/clc_client:latest'

  input:
    tuple val(dna_read1), val(dna_read2), val(rna_read1), val(rna_read2)
    // windows path, where the windows-hosted clc server expects to import files from
    // must be the full mounted path if its on a network drive and properly slash escaped
    // \\host.domain\MOUNTPOINT\SHARENAME\IMPORT\...
    // is therefore prefixed wth clc://serverfile
    val(clc_import_dir)
    // windows path, where the windows-hosted clc server expects to export files to
    // same rules as for clc_import_dir
    // can be the same as clc_import_dir
    val(clc_export_dir)
    // clc path where clc should write workflow outputs (clc storage, e.g. clc://server/...)
    val(clc_destdir)
    // string that specifies the clc workflow to call as required by the clc command line tools
    val(clc_workflow_name)
    // unix path, where the network share is mounted, to copy the results from there
    val(nas_import_dir)
    val(nas_export_dir)
    // workflow run id
    val(workflow_run_id)

  output:
    stdout emit: output
    val(vcf_path), emit: vcf_path
    val(csv_path), emit: csv_path

    script:

    def run_id_clc = workflow_run_id.replace(":","_")
    def workflow_run_destdir = "${clc_destdir + '/' + run_id_clc}"
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    def n3 = file(rna_read1).name
    def n4 = file(rna_read2).name
    //todo add checks for filen naming scheme, n1, ... n4, should have the same prefix
    // this is used below as the a directory to group the outputs related to the sample
    def samplename = samplename_from_filename(n1)

    // compose name of clc output files (this depends on the clc workflow settings)
    def output_filename = "${samplename}"

    vcf_path = "${nas_export_dir}/${samplename}/${output_filename}.vcf"
    csv_path = "${nas_export_dir}/${samplename}/${output_filename}.csv"

    // clc requires to specify output folders, iirc doesnt allow setting filename with the cli tools
    """
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A mkdir -t "${clc_destdir}" -n "${workflow_run_id}"
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A mkdir -t "${workflow_run_destdir}" -n "${samplename}"
    clcserver -S \$CLC_HOST -U \$CLC_USER -W \$CLC_PSW -A ${clc_workflow_name} \\
        --dna-reads-import-command ngs_import_illumina \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n1}\"  \\
        --dna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n2}\"  \\
        --rna-reads-import-algo-id ngs_import_illumina \\
        --rna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n3}\"  \\
        --rna-reads-select-files \"clc://serverfile/${clc_import_dir}/${n4}\"  \\
        --export-table-csv-csv-destination \"clc://serverfile/${clc_export_dir}/${samplename}\"  \\
        --export-table-csv-custom-file-name "${samplename}.csv\"  \\
        --export-vcf-vcf-destination \"clc://serverfile/${clc_export_dir}/${samplename}\"  \\
        --export-vcf-custom-file-name "${samplename}.vcf\"  \\
        -d "${workflow_run_destdir}/${samplename}"
    """

    stub:

    def destdir = clc_destdir
    def n1 = file(dna_read1).name
    def n2 = file(dna_read2).name
    def n3 = file(rna_read1).name
    def n4 = file(rna_read2).name
    def samplename = samplename_from_filename(n1)
    vcf_path = "${nas_import_dir}/${samplename}/${samplename}.vcf"
    csv_path = "${nas_import_dir}/${samplename}/${samplename}.csv"

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
  cache false

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

process mk_exportdir_dna_only {
  // copies files into an clc import/exportdir
  input:
    tuple val(dna_read1), val(dna_read2)
    val(nas_export_dir)
  output:
    tuple val(dna_read1), val(dna_read2)

  script:
  def n1 = file(dna_read1).name
  def samplename = samplename_from_filename(n1)
  def sample_exportdir = "${nas_export_dir}/${samplename}"

  """
  mkdir -p "${sample_exportdir}"
  """
}

process mk_exportdir_dna_rna {
  // copies files into an clc import/exportdir
  input:
    tuple val(dna_read1), val(dna_read2), val(rna_read1), val(rna_read2) 
    val(nas_export_dir)
  output:
    tuple val(dna_read1), val(dna_read2), val(rna_read1), val(rna_read2) 

  script:
  def n1 = file(dna_read1).name
  def samplename = samplename_from_filename(n1)
  def sample_exportdir = "${nas_export_dir}/${samplename}"

  """
  mkdir -p "${sample_exportdir}"
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



workflow pancancer_dna_only {
  take:
    samples
    clc_import_dir
    clc_export_dir
    clc_destdir
    clc_workflow_name
    nas_import_dir
    nas_export_dir
    workflow_run_id

  main:
    // check that samplenames match and are in the
    def check_row = { row -> 
        def sid = file(row.dna_read1).name.split('-')[0]
	def year = year_from_path(file(row.dna_read1).name)
	def expected_row = [
		"${sid}-${year}-DNA_1.fq.gz",
		"${sid}-${year}-DNA_2.fq.gz",
	]
	def actual_row = [
		file(row.dna_read1).name,
		file(row.dna_read2).name,
	]
	if (actual_row != expected_row){
		throw new Exception("row doesnt match naming scheme ${row}")
	}
    }

    def row_to_dict_dna_only = { row -> [dna_read1: row[0], dna_read2: row[1]]}

    samples.subscribe{ row -> check_row(row) }
    sample_files = samples.flatMap{ it -> ["${it.dna_read1}", "${it.dna_read2}"] }

    // use the files tuple for synchronizing the staging
    files = copyfiles(sample_files, nas_import_dir)
    staged_reads = files.map{it -> ["${samplename_from_path(it)}", it]}.groupTuple(by: 0, size: 2, sort: true).map{it -> it[1]}

    staged_reads.subscribe{ row -> check_row(row_to_dict_dna_only(row)) }
    staged_reads_w_export_dir = mk_exportdir_dna_only(staged_reads, nas_export_dir)

    out = clc_pancancer_dna_only(staged_reads_w_export_dir, 
    	clc_import_dir, clc_export_dir,
    	clc_destdir, clc_workflow_name,
	nas_import_dir, nas_export_dir,
	workflow_run_id)

    fils = copyback(out.vcf_path.mix(out.csv_path))

    fils.branch{
	vcf: it.getExtension() == 'vcf'
	csv: it.getExtension() == 'csv'
    }.set{ res }
  

    emit:
      vcf = res.vcf
      csv = res.csv

}

workflow pancancer_dna_rna {
  take:
    samples
    clc_import_dir
    clc_export_dir
    clc_destdir
    clc_workflow_name
    nas_import_dir
    nas_export_dir
    workflow_run_id

  main:
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

    def row_to_dict_dna_rna = { row ->
        [ dna_read1: row[0], dna_read2: row[1],
          rna_read1: row[2], rna_read2: row[3]
        ]
    }

    samples.subscribe{ row -> check_row(row) }

    sample_files = samples.flatMap{ it -> [
    	"${it.dna_read1}", "${it.dna_read2}", 
	    "${it.rna_read1}", "${it.rna_read2}"]
	}


    // use the files tuple for synchronizing the staging
    files = copyfiles(sample_files, nas_import_dir)
    staged_reads = files.map{it -> ["${samplename_from_path(it)}", it]}.groupTuple(by: 0, size: 4, sort: true).map{it -> it[1]}

    staged_reads.subscribe{ row -> check_row(row_to_dict_dna_rna(row)) }
    staged_reads_w_export_dir = mk_exportdir_dna_rna(staged_reads, nas_export_dir)

    out = clc_pancancer_dna_rna(staged_reads_w_export_dir, 
    	clc_import_dir, clc_export_dir,
    	clc_destdir, clc_workflow_name,
	nas_import_dir, nas_export_dir,
	workflow_run_id)

    fils = copyback(out.vcf_path.mix(out.csv_path))

    fils.branch{
	vcf: it.getExtension() == 'vcf'
	csv: it.getExtension() == 'csv'
    }.set{ res }
  

    emit:
      vcf = res.vcf
      csv = res.csv

}


workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  samples = Channel.fromPath(args.samplesheet, checkIfExists: true, type: 'file').splitCsv(header: true)

  pancancer_dna_only(samples,
  	args.clc_import_dir, 
	args.clc_export_dir,
	args.clc_destdir,
	args.clc_workflow_name,
	args.nas_import_dir,
	args.nas_export_dir
	)

  
  
}
