nextflow.enable.dsl=2

process bcfilter_and_annotsv {
   
    tag "${sample_id}"
    debug true

    //conda "bioconda::bcftools=1.21"
    
    publishDir(path: "${outdir}/${sample_id}", mode: "copy")

    input:
        tuple val(sample_id), path(sv_vcf)

    output:
        tuple val(sample_id), path("${sv_vcf.getSimpleName()}_filtered.vcf"), emit: filtered_sv_vcf
        tuple val(sample_id), path("${sv_vcf.getSimpleName()}_filtered_annotsv.tsv"), emit: sv_tsv
        tuple val(sample_id), path("${sv_vcf.getSimpleName()}_annotsv.log"), emit: annotsv_log 
        
    script:
    """   
    bcftools view -f PASS --output-type v ${sv_vcf} -o ${sv_vcf.getSimpleName()}_filtered.vcf     
    /projects/wgs_pilot/sv/AnnotSV/bin/AnnotSV -SVinputFile ${sv_vcf.getSimpleName()}_filtered.vcf \\
                                               -outputFile ${sv_vcf.getSimpleName()}_filtered_annotsv.tsv \\
                                               -genomeBuild GRCh38 \\
                                               -outputDir dir_${sv_vcf.getSimpleName()} > ${sv_vcf.getSimpleName()}_annotsv.log 2>&1

    mv dir_${sv_vcf.getSimpleName()}/${sv_vcf.getSimpleName()}_filtered_annotsv.tsv \\
        ${sv_vcf.getSimpleName()}_filtered_annotsv.tsv
    """
}

process ensembl_vep {

    tag "${sample_id}"
    debug true

    conda "bioconda::bcftools=1.21 bioconda::ensembl-vep=113.3"

    publishDir(path: "${outdir}/${sample_id}", mode: "copy")

    input:
       tuple val(sample_id), path(filtered_sv_vcf)

    output:
       tuple val(sample_id), path("${filtered_sv_vcf.getSimpleName()}_vep.ann.tab"), emit: vep_tab
       tuple val(sample_id), path("${filtered_sv_vcf.getSimpleName()}_vep.log"), emit: vep_log

    script:
    """
    vep \\
        -i ${filtered_sv_vcf} \\
        --everything \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        --format vcf \\
        -o ${filtered_sv_vcf.getSimpleName()}_vep.ann.tab \\
        --no_stats \\
        --cache \\
        --refseq \\
        --offline \\
        --dir_cache /projects/reference/vep_cache/113/indexed_refseq/ \\
        --fasta /projects/reference/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \\
        --cache_version 113 \\
        --terms 'SO' \\
        --distance 5000 \\
        --max_sv_size  1000000000 \\
        --fork 6 \\
        --pick_allele_gene \\
        --force_overwrite \\
        --overlaps \\
        --tab > ${filtered_sv_vcf.getSimpleName()}_vep.log > 2&>1
    """
}

process sv_format {

    tag "${sample_id}"
    debug true

    conda "conda-forge::python=3.9.18 conda-forge::pandas=2.1.0 conda-forge::openpyxl=3.1.2"

    publishDir(path: "${outdir}/sv_format_out", mode: "copy")

    input:
       tuple val(sample_id), path(sv_tsv), path(vep_tab)
       val(genelist)

    output:
       tuple val(sample_id), path("*.sv.csv"), emit: sv_csv
       //tuple val(sample_id), path("*.annoformat.log"), emit: annoformat_log      

    script:
    """
    grep -v "##" ${vep_tab} > ${vep_tab.getSimpleName()}.noheader.tab
    sv_annoformat.py \\
        --sv_table ${sv_tsv} \\
        --wgs_pilot_gene_list ${genelist} \\
        --vep_tab ${vep_tab.getSimpleName()}.noheader.tab \\
        --sample ${sample_id} \\
        --outfile ${sample_id}.sv.csv
   """
}

csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, file(row.sv_vcf)) }
outdir = params.outdir
genelist = Channel.value(params.genelist)

workflow prepare_sv {

    take:
       csv_ch
       genelist

   main:
      bcfilter_and_annotsv(csv_ch)
      ensembl_vep(bcfilter_and_annotsv.out.filtered_sv_vcf)
      sv_format_ch = bcfilter_and_annotsv.out.sv_tsv.join(ensembl_vep.out.vep_tab)
      sv_format(sv_format_ch,genelist)
}

workflow {
    prepare_sv(csv_ch,genelist)
}
