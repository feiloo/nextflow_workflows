include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"
include { clc_nextflow } from "$NEXTFLOW_MODULES/clc_nextflow"
include { CIOABCD_VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation"
include { sequence_alignment } from "$NEXTFLOW_MODULES/sequence_alignment"
include { pancancer_analyse } from "$NEXTFLOW_MODULES/nextpipe"

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
	clc_out = clc_nextflow(args.samplesheet, 
	  args.clc_import_dir, 
	  args.clc_export_dir,
	  args.clc_destdir,
	  args.clc_workflow_name,
	  args.nas_import_dir,
	  args.nas_export_dir,
	  args.workflow_run_id
	  )

	pancancer_analyse(
	  clc_out.vcf, clc_out.csv, args.vep_cache, args.vep_refgenome, 
	  args.transcriptlist, args.variantlist, args.outdir
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
	new_rows = Channel.of("sample,vcf").concat((fixed_vcfs.map{it -> "${it[0]},${it[1]}"})).collect()
	//.collectFile(name: "fixed_samplesheet.csv", newLine: true)
	new_samplesheet = write_samplesheet(new_rows).samplesheet
	fasta = rename_chromosomes_refgenome(args.refgenome)

  	CIOABCD_VARIANTINTERPRETATION(args, new_samplesheet, fasta)
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

  	CIOABCD_VARIANTINTERPRETATION(args, new_samplesheet, fasta)

	publish(pub, args.output_dir)
  }
  else {
    println "invalid value for parameter workflow_variation " + args.workflow_variation
    System.exit(1)
  }
}
