nextflow.enable.dsl=2

process bam_index {
   
    tag "${sample_id}"
    debug true

    conda "bioconda::samtools=1.21"

    input:
        tuple val(sample_id), val(female), path(normal_bam), path(tumor_bam)

    output:
        tuple val(sample_id), val(female), path("${normal_bam}"), path("${normal_bam}.bai"), path("${tumor_bam}"), path("${tumor_bam}.bai"), emit: bam_bai

    script:
    def args = task.ext.args ?: ""                                                                                                                                               
    def args2 = task.ext.args2 ?: "" 
    """
    samtools index -@ $args ${normal_bam} ${normal_bam}.bai
    samtools index -@ $args2 ${tumor_bam} ${tumor_bam}.bai
    """
}

process bam2seqz {
    
    tag "${sample_id}"
    debug true
            
    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7' 

    //publishDir(path: "${outdir}/${sample_id}/preprocessed_files", mode: "copy")

    input:
        tuple val(sample_id), val(female), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
        val(ref_fasta_gz)
        val(gc_wig_file)
        
    output:
        tuple val(sample_id), val(female), path("${sample_id}*.seqz.gz"), emit: seqzfile
        tuple val(sample_id), val(female), path("${sample_id}*.seqz.gz.tbi"), emit: tbi_seqzfile
            
    script:
    def args = task.ext.args ?: ""     
    def args2 = task.ext.args2 ?: ""
    """
    sequenza-utils bam2seqz \\
                   --normal ${normal_bam} \\
                   --tumor ${tumor_bam} \\
                   --fasta ${ref_fasta_gz} \\
                   -C $args \\
                   --parallel $args2 \\
                   -gc ${gc_wig_file} \\
                   --output ${sample_id}.seqz.gz
    """
}

process merge_seqz_files {

    tag "${sample_id}"
    debug true
     
    //publishDir(path: "${outdir}/${sample_id}/merged_seqz_file", mode: "copy")

    input:
        tuple val(sample_id), val(female), path(seqzfile)    

    output:
        tuple val(sample_id), val(female), path("header_${sample_id}.txt"),     emit: header_seqzfile
        tuple val(sample_id), val(female), path("merged_${sample_id}.seqz"),    emit: merged_seqzfile
        
     script:
     """
     for file in `ls -v ${seqzfile}` 
     do 
         zcat \${file} | awk 'NR==1 {print; exit}' > header_${sample_id}.txt
         zcat \${file} | awk 'NR>1 {print}' >> merged_${sample_id}.seqz
     done 
     """
}

process formatting_merged_seqz_file {

    tag "${sample_id}"
    debug true

    //publishDir(path: "${outdir}/${sample_id}/merged_seqz_file", mode: "copy")

    input:
        tuple val(sample_id), val(female), path(header_seqzfile)
        tuple val(sample_id), val(female), path(merged_seqzfile)
        
    output:
        tuple val(sample_id), val(female), path("all_${sample_id}.seqz.gz"), path("all_${sample_id}.seqz.gz.tbi")

     script:
     def args = task.ext.args ?: ""   
     """
     cat ${header_seqzfile} ${merged_seqzfile} | bgzip -@ $args > all_${sample_id}.seqz.gz  
     tabix -f -s 1 -b 2 -e 2 -S 1 all_${sample_id}.seqz.gz
     """

}

process seqz_binning {
    
    tag "${sample_id}"
    debug true

    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7'

    //publishDir(path: "${outdir}/${sample_id}/preprocessed_files", mode: "copy")
    
    input:
        tuple val(sample_id), val(female), path(all_seqzfile), path(all_tbi_seqzfile)

    output:
        tuple val(sample_id), val(female), path("${sample_id}_bin.seqz.gz"), emit: smallbinfile

    script:
    def args = task.ext.args ?: ""
    """
    sequenza-utils seqz_binning \\
              --seqz ${all_seqzfile} \\
              --window $args \\
              -o ${sample_id}_bin.seqz.gz
    """
}

process r_sequenza { 
     
    tag "${sample_id}"
    debug true

    //conda "bioconda::r-sequenza=3.0.0=r40h3342da4_3"
    
    publishDir(path: "${outdir}/${sample_id}/sequenza_results", mode: "copy")
    
    input:
       tuple val(sample_id), val(female), path(smallbinfile)
  
    output:
       tuple val(sample_id), path("${sample_id}_alternative_fit.pdf"),        emit: alternative_fit 
       tuple val(sample_id), path("${sample_id}_alternative_solutions.txt"),  emit: alternative_solutions
       tuple val(sample_id), path("${sample_id}_chromosome_depths.pdf"),      emit: chromosome_depths
       tuple val(sample_id), path("${sample_id}_chromosome_view.pdf"),        emit: chromosome_view
       tuple val(sample_id), path("${sample_id}_CN_bars.pdf"),                emit: CN_bars
       tuple val(sample_id), path("${sample_id}_confints_CP.txt"),            emit: confints_CP
       tuple val(sample_id), path("${sample_id}_CP_contours.pdf"),            emit: CP_contours
       tuple val(sample_id), path("${sample_id}_gc_plots.pdf"),               emit: gc_plots
       tuple val(sample_id), path("${sample_id}_genome_view.pdf"),            emit: genome_view
       tuple val(sample_id), path("${sample_id}_model_fit.pdf"),              emit: model_fit
       tuple val(sample_id), path("${sample_id}_mutations.txt"),              emit: mutations
       tuple val(sample_id), path("${sample_id}_segments.txt"),               emit: segments
       tuple val(sample_id), path("${sample_id}_sequenza_cp_table.RData"),    emit: sequenza_cp_table
       tuple val(sample_id), path("${sample_id}_sequenza_extract.RData"),     emit: sequenza_extract
       tuple val(sample_id), path("${sample_id}_sequenza_log.txt"),           emit: sequenza_log
       tuple val(sample_id), path("${sample_id}_cellularity_ploidy.txt"),     emit: cellularity_ploidy

    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    data.file <- "${smallbinfile}"
    seqzdata <- sequenza.extract(data.file, assembly="hg38")
    CP <- sequenza.fit(seqzdata, female="${female}")
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${sample_id}", out.dir = ".")
    
    cint <- get.ci(CP)
    cellularity <- cint[["max.cellularity"]]
    ploidy <- cint[["max.ploidy"]]
    data_CP <- data.frame(cellularity = cellularity, ploidy = ploidy)
    write.table(data_CP, file="${sample_id}_cellularity_ploidy.txt", sep="\t", row.names=FALSE)
    """
}     

process preprare_scar_hrd_input {

    tag "${sample_id}"
    debug true

    input:
      tuple val(sample_id), path(segments)
      tuple val(sample_id), path(cellularity_ploidy)
      val(scar_hrd_header)

    output:
       tuple val(sample_id), path("${sample_id}_scarHRD_input.txt"), emit: scarhrd_input

    script:
    """
    awk 'BEGIN {FS=OFS="\t"} {print \$1,\$2,\$3,\$10,\$11,\$12}' ${segments} > tmp1_${sample_id}.txt
    sed -i "s/\$/\t${sample_id}/" tmp1_${sample_id}.txt
    ploidyValue=`awk 'NR>1{print \$2}' ${cellularity_ploidy}`
    sed -i "s/\$/\t\$ploidyValue/" tmp1_${sample_id}.txt
    awk 'BEGIN {FS=OFS="\t"} {print \$7,\$1,\$2,\$3,\$4,\$5,\$6,\$8}' tmp1_${sample_id}.txt > tmp2_${sample_id}.txt
    awk 'NR>1 {print}' tmp2_${sample_id}.txt > tmp3_${sample_id}.txt
    cat ${scar_hrd_header} tmp3_${sample_id}.txt > ${sample_id}_scarHRD_input.txt 
    """
}

process scar_hrd {
    
    tag "${sample_id}"
    debug true

    publishDir(path: "${outdir}/${sample_id}/scarHRD_results", mode: "copy")

    input:
      tuple val(sample_id), path(scarhrd_input)

    output:
      tuple val(sample_id), path("${sample_id}_HRDresults.txt"), emit: scarHRDfile

     script:
     def args = task.ext.args ?: ""
     """
     #!/usr/bin/env Rscript
     library(scarHRD)
     scarHRD_input <- "${scarhrd_input}"
     scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=$args)
     """
}


csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, row.female, file(row.normal_bam), file(row.tumor_bam)) }
ref_fasta_gz_ch = Channel.value(params.ref_fasta_gz)
gc_wig_file_ch = Channel.value(params.gc_wig_file)
outdir = params.outdir
scar_hrd_header_ch = Channel.value(params.scar_hrd_header)

workflow sequenza {

    take:
       csv_ch
       ref_fasta_gz_ch
       gc_wig_file_ch
       scar_hrd_header_ch

       
   main:
       bam_index(csv_ch)
       bam2seqz(bam_index.out.bam_bai,ref_fasta_gz_ch,gc_wig_file_ch)
       merge_seqz_files(bam2seqz.out.seqzfile)
       formatting_merged_seqz_file(merge_seqz_files.out.header_seqzfile,merge_seqz_files.out.merged_seqzfile)
       seqz_binning(formatting_merged_seqz_file.out)
       r_sequenza(seqz_binning.out.smallbinfile)
       preprare_scar_hrd_input(r_sequenza.out.segments,r_sequenza.out.cellularity_ploidy,scar_hrd_header_ch)
       scar_hrd(preprare_scar_hrd_input.out.scarhrd_input)
}

workflow {
    sequenza(csv_ch,ref_fasta_gz_ch,gc_wig_file_ch,scar_hrd_header_ch)
}
