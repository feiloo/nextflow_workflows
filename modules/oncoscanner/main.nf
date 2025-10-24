include { arriba_nextflow } from "$NEXTFLOW_MODULES/arriba_nextflow"

if(params.workflow_variation == 'clc'){
include { pancancer_dna_only; pancancer_dna_rna } from "$NEXTFLOW_MODULES/clc_nextflow"
include { pancancer_analyse } from "$NEXTFLOW_MODULES/nextpipe"
}

include { VARIANTINTERPRETATION } from "$NEXTFLOW_MODULES/variantinterpretation/workflows/variantinterpretation.nf"
include { sequence_alignment } from "$NEXTFLOW_MODULES/sequence_alignment"
include { analyse_biomarkers } from "$NEXTFLOW_MODULES/biomarker"

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

process bgzip_vcf {
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"
    cpus 4

    input:
    path(vcf)

    output:
    path("out/${vcf.getSimpleName()}.vcf.gz"), emit: vcf
 
    script:
    """
    mkdir out
    bgzip -@ ${task.cpus} "${vcf}" -c >"out/${vcf.getSimpleName()}.vcf.gz"
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
        samples = Channel.fromPath(args.samplesheet, checkIfExists: true, type: 'file').splitCsv(header: true)

	untared_vep = untar_file(args.vep_cache)

        if(args.containsKey('has_rna') && args.has_rna == true) {
          clc_out = pancancer_dna_rna(
            samples,
            args.clc_import_dir, 
            args.clc_export_dir,
            args.clc_destdir,
	    args.clc_dna_rna_workflow_name,
            args.nas_import_dir,
            args.nas_export_dir,
            args.workflow_run_id
            )
        } else {
          clc_out = pancancer_dna_only(
            samples,
            args.clc_import_dir, 
            args.clc_export_dir,
            args.clc_destdir,
	    args.clc_dna_only_workflow_name,
            args.nas_import_dir,
            args.nas_export_dir,
            args.workflow_run_id
            )
        }

	pancancer_analyse(
	  clc_out.vcf, clc_out.csv, untared_vep, args.vep_refgenome, 
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
	  args.library_type,
	  )
	
  }
  else if (args.workflow_variation == 'align_interpret'){
	seq_output = sequence_alignment(args)
	// preprocessing output
	pub = seq_output.integrity_check.mix(seq_output.samplesheet)
	pub = seq_output.samplesheet.mix(seq_output.hash_db)
	// fastp report currently broken, skipping it
	pub = pub.mix(seq_output.fastp_report_h).mix(seq_output.fastp_report_j)
	pub = pub.mix(seq_output.vcf).mix(seq_output.bam_coverage).mix(seq_output.bam_stats)

	//pub = pub.mix(seq_output.bam_pairs_w_idx.flatten())//.mix(seq_output.vcf).mix(seq_output.bam_coverage).mix(seq_output.bam_stats)

	biomarkers = analyse_biomarkers(seq_output.bam_pairs_w_idx, args.refgenome, seq_output.refgenome_index, seq_output.refgenome_dict, args)
	msi_csv = biomarkers.msi_csv
	hs_metrics = biomarkers.hs_metrics
	// matched_bams_w_idx = biomarkers.matched_preproc_bams

	pub = pub.mix(msi_csv).mix(hs_metrics) //.mix(matched_bams_w_idx.flatten())
	pub = pub.mix(seq_output.bam_w_idx.flatten())

	matched_vcf = bgzip_vcf(seq_output.vcf).vcf.map{it -> [it.getSimpleName(), it]}

	chr_fixed_vcf = rename_chromosomes_vcf(matched_vcf).vcf
	fixed_vcfs = rename_clcad_to_ad(chr_fixed_vcf).vcf
	new_rows = Channel.of("sample,vcf").concat(fixed_vcfs.map{it -> "${it[0]},${it[1]}"}).collect()
	//.collectFile(name: "fixed_samplesheet.csv", newLine: true)
	new_samplesheet = write_samplesheet(new_rows).samplesheet
	fasta = rename_chromosomes_refgenome(args.refgenome)

	extra_files = []
	annotation_colinfo = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/annotation_colinfo.tsv")
	datavzrd_config = Channel.fromPath("$NEXTFLOW_MODULES/variantinterpretation/assets/datavzrd_config_template.yaml")

	new_csv_channel = fixed_vcfs.map { it -> [[id:it[0]], it[1]] }

	untared_vep = untar_file(args.vep_cache)

  	interpretation = VARIANTINTERPRETATION(
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
	  args.library_type,
	  ).ukb_results
	
	pub = pub.mix(interpretation)

	def skip = args.skip_publishing?.toString()?.toLowerCase() ?: 'false'
	if (skip == 'false' || skip == '') {
		publish(pub, args.output_dir)
	} else if (skip == 'true'){
		// do not publish anything
		println "skipping publishing"
	} else {
	    println "invalid value for parameter skip_publishing " + args.workflow_variation
	    System.exit(1)
	}
  }
  else {
    println "invalid value for parameter workflow_variation " + args.workflow_variation
    System.exit(1)
  }
}

workflow.onComplete {
    println "Pipeline for samplesheet $params.samplesheet finished with outputs at: $params.output_dir"
    println "Pipeline used workdir $params.workdir"
}
