nextflow.enable.dsl=2

include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"
include { clc_nextflow } from "$NEXTFLOW_MODULES/clc_nextflow"
include { CIOABCD_VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation"
include { sequence_alignment } from "$NEXTFLOW_MODULES/sequence_alignment"
include { msi_annotate } from "$NEXTFLOW_MODULES/biomarker"
include { gatk_collect_hs_metrics; gatk_bed_to_intervallist; gatk_createsequencedictionary } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"
include { index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
//include { SAREK } from "$NEXTFLOW_MODULES/sarek_wrapper"

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

process get_sha256sum {
    input:
    tuple val(sample), path(vcf), val(samplesheet_sha256sum)

    output:
    tuple val(sample), path(vcf), val(samplesheet_sha256sum), stdout
    
    script:
    """
    sha256sum "${vcf}" | cut -d ' ' -f 1
    """

}

process save_samplesheet {
    publishDir "${args.output_dir}", mode: 'copy', overwrite: false

    input:
      path(samplesheet)
    output:
      path(samplesheet)
    script:
      println "saved samplesheet to output_dir: ${args.output_dir}"
}


workflow {
  //def args = params
  def args = [:]
  for (param in params) { args[param.key] = param.value }

  if(args.workflow_variation == 'sequence_alignment'){
  	sequence_alignment(args)
  }
  else if(args.workflow_variation == 'arriba'){
  	arriba_nextflow(args)
  }
  else if(args.workflow_variation == 'clc'){
  	clc_nextflow(args)
  }
  else if(args.workflow_variation == 'sarek'){
  	SAREK(args)
  }
  else if(args.workflow_variation == 'msisensorpro'){
	samplesheet = args.samplesheet
	header = ['sample', 'tumor_bam', 'tumor_sha256sum', 'normal_bam', 'normal_sha256sum']
    	rows = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: header, skip: 1)
	bams = rows.map{it -> [it.tumor_bam, it.normal_bam]}
	msi_annotate(bams, args.refgenome)

	tumor_bams = bams.map{it -> [it[0]]}

	refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict
	refgenome_index = index_fasta(args.refgenome).fasta_index

	targets_list = gatk_bed_to_intervallist(args.targets_bed, args.refgenome, refgenome_dict).targets_list
	gatk_collect_hs_metrics(tumor_bams, targets_list, args.refgenome, refgenome_index)

  }
  else if(args.workflow_variation == 'variantinterpretation'){
	samplesheet = args.samplesheet
	header = ['sample', 'vcf', 'vcf_sha256sum']
    	rows = Channel.fromPath(samplesheet, checkIfExists: true, type: 'file').splitCsv(header: header, skip: 1)

	rows_w_shasums = get_sha256sum(rows)

	save_samplesheet(samplesheet)

	correct_rows = rows_w_shasums.map{it -> 
		if("${it[2]}".trim() != "${it[3]}".trim()) {
		  throw new Exception("samplesheet hash doesnt match actual filehash ${it}")
		  }
		else {
		  [it[0], it[1]]
		  }
		}

	//alternatively add a .filter op here, if we wanna fail late

	chr_fixed_vcf = rename_chromosomes_vcf(correct_rows).vcf
	fixed_vcfs = rename_clcad_to_ad(chr_fixed_vcf).vcf
	new_rows = Channel.of("sample,vcf").concat((fixed_vcfs.map{it -> "${it[0]},${it[1]}"})).collect()
	//.collectFile(name: "fixed_samplesheet.csv", newLine: true)
	new_samplesheet = write_samplesheet(new_rows).samplesheet
	fasta = rename_chromosomes_refgenome(args.refgenome)

  	CIOABCD_VARIANTINTERPRETATION(args, new_samplesheet, fasta)
  }
  else {
    println "invalid value for parameter workflow_variation " + args.workflow_variation
    System.exit(1)
  }
}
