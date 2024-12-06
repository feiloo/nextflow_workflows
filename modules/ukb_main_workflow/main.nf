include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"
include { clc_nextflow } from "$NEXTFLOW_MODULES/clc_nextflow"
include { VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation/workflows/variantinterpretation.nf"
include { sequence_alignment } from "$NEXTFLOW_MODULES/sequence_alignment"

include { publish } from "$NEXTFLOW_MODULES/sequence_alignment/utils.nf"

process bcftools_get_tumor {
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"

    input:
    path(vcf)

    output:
    path("out/${vcf.getSimpleName()}.vcf.gz"), emit: vcf
 
    script:
    """
    mkdir out
    bcftools view --samples "${vcf.getSimpleName().split('_')[0]}_T_1" "${vcf}" | bgzip -c - >"out/${vcf.getSimpleName()}.vcf.gz"
    """
}

process rename_clcad_to_ad {
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("out/${vcf}"), emit: vcf
 
    script:
    """
    mkdir out
    bgzip -cd "${vcf}" | sed -e 's/CLCAD2/AD/g' - | bgzip -c > "out/${vcf}"
    """
}

process rename_chromosomes_vcf {
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"

    fair true

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("out/${vcf}"), emit: vcf

    script:
    def chr_map = (1..22).collect{"${it} chr${it}"}
    chr_map += ["X chrX", "Y chrY"]
    chr_map = chr_map.join('\n')
    """
    echo "${chr_map}" > chr_map.txt
    mkdir out
    bcftools annotate --rename-chrs chr_map.txt --output out/${vcf} ${vcf}
    """
}

process rename_chromosomes_refgenome {
    //
    input:
    path(refgenome)

    output:
    path("out/${refgenome}")

    script:
    """
    mkdir mid
    mkdir out
    sed 's/^>\\([1-9XY]\\)/>chr\\1/g' "${refgenome}" > "mid/${refgenome}"
    dos2unix -F -n "mid/${refgenome}" "out/${refgenome}"
    
    """
}

process write_samplesheet {
    input:
    val(samplesheet_rows)

    output:
    path("fixed_samplesheet.csv"), emit: samplesheet

    script:
    def samplesheet_text = samplesheet_rows.join('\n')
    """
    echo "${samplesheet_text}" > fixed_samplesheet.csv
    """
}

process untar_file {
    input:
    path(filen)

    output:
    path(untared_cache)

    script:
    """
    mkdir untared_cache
    tar -xzf ${filen} --directory untared_cache
    """
}


workflow {
  //def args = params
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  if(args.workflow_variation == 'sequence_alignment'){
	output = sequence_alignment(args)
	pub = output.bam.mix(output.vcf).mix(output.bam_coverage).mix(output.bam_stats)
	publish(pub, args.output_dir)
  }
  else if(args.workflow_variation == 'arriba'){
  	arriba_nextflow(args)
  }
  else if(args.workflow_variation == 'clc'){
  	//clc_nextflow(args)

	  clc_nextflow(args.samplesheet, 
		args.clc_import_dir, 
		args.clc_export_dir,
		args.clc_destdir,
		args.clc_workflow_name,
		args.nas_import_dir,
		args.nas_export_dir
		)
  }
  else if(args.workflow_variation == 'sarek'){
  	SAREK(args)
  }
  else if(args.workflow_variation == 'variantinterpretation'){
	samplesheet = args.samplesheet
	header = ['sample', 'vcf']
    	csv_channel = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: header, skip: 1)
	chr_fixed_vcf = rename_chromosomes_vcf(csv_channel).vcf
	fixed_vcfs = rename_clcad_to_ad(chr_fixed_vcf).vcf
	//new_rows = Channel.of("sample,vcf").concat((fixed_vcfs.map{it -> "${it[0]},${it[1]}"})).collect()

	//.collectFile(name: "fixed_samplesheet.csv", newLine: true)

	//new_samplesheet = write_samplesheet(new_rows).samplesheet

	fasta = rename_chromosomes_refgenome(args.refgenome)

	extra_files = []
	annotation_colinfo = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/annotation_colinfo.tsv")
	datavzrd_config = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/datavzrd_config_template.yaml")

	new_csv_channel = fixed_vcfs.map { it -> [[id:it[0]], it[1]] }

	untared_vep = untar_file(args.vep_cache)

  	VARIANTINTERPRETATION(
	  new_csv_channel, //new_samplesheet,
	  fasta,
	  untared_vep,
	  args.vep_cache_version,
	  args.vep_genome,
	  args.vep_species,
	  extra_files,
	  args.annotation_fields,
	  args.transcriptlist,
	  datavzrd_config,
	  annotation_colinfo,
	  args.bedfile,
	  args.custom_filters,
	  )
	
  }
  else if (args.workflow_variation == 'align_interpret'){
	output = sequence_alignment(args)
	pub = output.bam.mix(output.vcf).mix(output.bam_coverage).mix(output.bam_stats)

	tumor_vcf = bcftools_get_tumor(output.vcf).vcf.map{it -> [it.getSimpleName(), it]}

	chr_fixed_vcf = rename_chromosomes_vcf(tumor_vcf).vcf
	fixed_vcfs = rename_clcad_to_ad(chr_fixed_vcf).vcf
	new_rows = Channel.of("sample,vcf").concat((fixed_vcfs.map{it -> "${it[0]},${it[1]}"})).collect()
	//.collectFile(name: "fixed_samplesheet.csv", newLine: true)
	new_samplesheet = write_samplesheet(new_rows).samplesheet
	fasta = rename_chromosomes_refgenome(args.refgenome)

	extra_files = []
	annotation_colinfo = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/annotation_colinfo.tsv")
	datavzrd_config = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/datavzrd_config_template.yaml")

	new_csv_channel = fixed_vcfs.map { it -> [[id:it[0]], it[1]] }

	untared_vep = untar_file(args.vep_cache)

  	VARIANTINTERPRETATION(
	  new_csv_channel, //new_samplesheet,
	  fasta,
	  untared_vep,
	  args.vep_cache_version,
	  args.vep_genome,
	  args.vep_species,
	  extra_files,
	  args.annotation_fields,
	  args.transcriptlist,
	  datavzrd_config,
	  annotation_colinfo,
	  args.bedfile,
	  args.custom_filters,
	  )
	

	publish(pub, args.output_dir)
  }
  else {
    println "invalid value for parameter workflow_variation " + args.workflow_variation
    System.exit(1)
  }
}
