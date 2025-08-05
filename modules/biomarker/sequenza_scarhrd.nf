nextflow.enable.dsl=2

process prepare_bam_for_sequenza {
            
    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6 bioconda::r-sequenza=3.0.0 conda-forge::r-devtools=2.4.5"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7' 

    //publishDir(path: "${outdir}/${tumor_bam.getSimpleName()}/preprocessed_files", mode: "copy")
    cpus 24
    errorStrategy 'retry'
    memory '250 GB'

    input:
        tuple path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
        path(ref_fasta_gz)
        path(gc_wig_file)
        
    output:
        //path("${normal_bam.getSimpleName()}*.seqz.gz"), emit: seqzfile
        //path("${normal_bam.getSimpleName()}*.seqz.gz.tbi"), emit: tbi_seqzfile
        //path("header_${normal_bam.getSimpleName()}.txt"),     emit: header_seqzfile
        //path("merged_${normal_bam.getSimpleName()}.seqz"),    emit: merged_seqzfile
        //path("all_${sample_id}.seqz.gz"), path("all_${sample_id}.seqz.gz.tbi")
        //tuple val(sample_id), val(female), path(all_seqzfile), path(all_tbi_seqzfile)

        //tuple val(sample_id), path("${sample_id}_bin.seqz.gz"), emit: smallbinfile
        tuple val(sample_id), path("${sample_id}*.seqz.gz"), emit: seqzfile
        tuple val(sample_id), path("${sample_id}*.seqz.gz.tbi"), emit: tbi_seqzfile
    
    script:
    def args = task.ext.args ?: "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"     
    def args2 = task.ext.args2 ?: ""
    sample_id = "${normal_bam.getSimpleName()}"
    """
    sequenza-utils bam2seqz \\
                   --normal ${normal_bam} \\
                   --tumor ${tumor_bam} \\
                   --fasta ${ref_fasta_gz} \\
                   -C $args \\
                   --parallel ${task.cpus} \\
                   -gc ${gc_wig_file} \\
                   --output ${sample_id}.seqz.gz
    """
}

process merge_seqz_files {

    tag "${sample_id}"
    debug true
     
    //publishDir(path: "${outdir}/${sample_id}/merged_seqz_file", mode: "copy")

    input:
        tuple val(sample_id), path(seqzfile)    

    output:
        tuple val(sample_id), path("header_${sample_id}.txt"),     emit: header_seqzfile
        tuple val(sample_id), path("merged_${sample_id}.seqz"),    emit: merged_seqzfile
        
     script:
     def header_seqz="chromosome      position        base.ref        depth.normal    depth.tumor     depth.ratio     Af      Bf      zygosity.normal GC.percent      good.reads      AB.normal   AB.tumor tumor.strand"
     """
     echo $header_seqz > header_${sample_id}.txt
     for file in `ls -v ${seqzfile}` 
     do 
         zcat \${file} | awk 'NR>1 {print}' >> merged_${sample_id}.seqz
     done 
     """
}

process formatting_merged_seqz_file {

    tag "${sample_id}"
    debug true

    //publishDir(path: "${outdir}/${sample_id}/merged_seqz_file", mode: "copy")

    input:
        tuple val(sample_id), path(header_seqzfile)
        tuple val(sample_id), path(merged_seqzfile)
        
    output:
        tuple val(sample_id), path("all_${sample_id}.seqz.gz"), path("all_${sample_id}.seqz.gz.tbi")

     script:
     def args = task.ext.args ?: ""   
     """
     cat ${header_seqzfile} ${merged_seqzfile} | bgzip -@ 24 > all_${sample_id}.seqz.gz  
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
        tuple val(sample_id), path(all_seqzfile), path(all_tbi_seqzfile)

    output:
        tuple val(sample_id), path("${sample_id}_bin.seqz.gz"), emit: smallbinfile

    script:
    def args = task.ext.args ?: ""
    """
    sequenza-utils seqz_binning \\
              --seqz ${all_seqzfile} \\
              --window 50 \\
              -o ${sample_id}_bin.seqz.gz
    """
}

process r_sequenza { 

    tag "${sample_id}"
    debug true
    
    conda "/home/pbasitta/miniconda3/envs/r_seq"
   
    
    //publishDir(path: "${outdir}/${smallbinfile.getSimpleName()}/sequenza_results", mode: "copy")
    errorStrategy 'retry'
    
    input:
       tuple val(sample_id), path(smallbinfile), val(female)
  
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
       tuple val(sample_id), path("${sample_id}_sequenza_cp_table.RData"),    emit: sequenza_cp_table
       tuple val(sample_id), path("${sample_id}_sequenza_extract.RData"),     emit: sequenza_extract
       tuple val(sample_id), path("${sample_id}_sequenza_log.txt"),           emit: sequenza_log
       tuple val(sample_id), path("${sample_id}_segments.txt"),               emit: segments
       tuple val(sample_id), path("${sample_id}_cellularity_ploidy.txt"),     emit: cellularity_ploidy

       //tuple val("${smallbinfile.getSimpleName()}"), path("${smallbinfile.getSimpleName()}_segments.txt"), path("${smallbinfile.getSimpleName()}_cellularity_ploidy.txt"),     emit: segments_ploidy
       //path("${smallbinfile.getSimpleName()}_sequenza_results.tar.gz}"), emit: sequenza_results

    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    data.file <- "${smallbinfile}"
    seqzdata <- sequenza.extract(data.file, assembly="hg38")
    
    sex <- read.table(file="${female}", header=FALSE)
    if (length(sex) == 1) {
      if (sex == "female") {
	female <- TRUE
      } else if ( sex == "male") {
	female <- FALSE
      } else {
	stop("File contains unexpected content.")
      }
    } else {
      stop("File does not contain exactly one line.")
    }
    
    CP <- sequenza.fit(seqzdata, female=female)
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${sample_id}", out.dir = ".")
    
    cint <- get.ci(CP)
    cellularity <- cint[["max.cellularity"]]
    ploidy <- cint[["max.ploidy"]]
    data_CP <- data.frame(cellularity = cellularity, ploidy = ploidy)
    write.table(data_CP, file="${sample_id}_cellularity_ploidy.txt", sep="\t", row.names=FALSE)
    """

    //"""
    //mkdir sequenza_results
    //cp *.pdf sequenza_results
    //cp *.txt sequenza_results
    //cp *.RData sequenza_results
    //tar czf ${smallbinfile.getSimpleName()}_sequenza_results.tar.gz
    //"""
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
    //def scar_hrd_header = "SampleID	Chromosome	Start_position	End_position	total_cn	A_cn	B_cn	ploidy"
    """
    awk 'BEGIN {FS=OFS="\t"} {print \$1,\$2,\$3,\$10,\$11,\$12}' ${segments} > tmp1_${sample_id}.txt
    sed -i "s/\$/\t${sample_id}/" tmp1_${sample_id}.txt
    ploidyValue=`awk 'NR>1{print \$2}' ${cellularity_ploidy}`
    sed -i "s/\$/\t\$ploidyValue/" tmp1_${sample_id}.txt
    awk 'BEGIN {FS=OFS="\t"} {print \$7,\$1,\$2,\$3,\$4,\$5,\$6,\$8}' tmp1_${sample_id}.txt > tmp2_${sample_id}.txt
    awk 'NR>1 {print}' tmp2_${sample_id}.txt > tmp3_${sample_id}.txt
    cat $scar_hrd_header tmp3_${sample_id}.txt > ${sample_id}_scarHRD_input.txt 
    """
}

process scar_hrd {
    
    tag "${sample_id}"
    debug true
    
    conda "/home/pbasitta/miniconda3/envs/r_seq"

    //publishDir(path: "${outdir}/${sample_id}/scarHRD_results", mode: "copy")

    input:
      tuple val(sample_id), path(scarhrd_input)

    output:
      tuple val(sample_id), path("${sample_id}_HRDresults.txt"), emit: scarHRDfile

     script:
     """
     #!/usr/bin/env Rscript
     library(scarHRD)
     scarHRD_input <- "${scarhrd_input}"
     scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=TRUE)
     """
}

process scar_hrd_dev {
    
    tag "${sample_id}"
    debug true
    errorStrategy 'retry'

    conda "/home/pbasitta/miniconda3/envs/r_seq"

    //publishDir(path: "${outdir}/${sample_id}/scarHRD_results", mode: "copy")

    input:
      //tuple val(sample_id), path(segments), path(cellularity_ploidy)
      tuple val(sample_id), path(segments)
      tuple val(sample_id), path(cellularity_ploidy)

    output:
      path("${segments.getSimpleName()}_HRDresults.txt"), emit: scar_hrd_results

    script:
    def args = task.ext.args ?: ""
    def scar_hrd_header = 'SampleID        Chromosome      Start_position  End_position    total_cn        A_cn    B_cn    ploidy'

    """
    awk 'BEGIN {FS=OFS="\t"} {print \$1,\$2,\$3,\$10,\$11,\$12}' ${segments} > tmp1_${segments.getSimpleName()}.txt
    sed -i "s/\$/\t${segments.getSimpleName()}/" tmp1_${segments.getSimpleName()}.txt
    ploidyValue=`awk 'NR>1{print \$2}' ${cellularity_ploidy}`
    sed -i "s/\$/\t\$ploidyValue/" tmp1_${segments.getSimpleName()}.txt
    awk 'BEGIN {FS=OFS="\t"} {print \$7,\$1,\$2,\$3,\$4,\$5,\$6,\$8}' tmp1_${segments.getSimpleName()}.txt > tmp2_${segments.getSimpleName()}.txt
    awk 'NR>1 {print}' tmp2_${segments.getSimpleName()}.txt > tmp3_${segments.getSimpleName()}.txt
    echo $scar_hrd_header > scar_hrd_header.txt
    cat scar_hrd_header.txt tmp3_${segments.getSimpleName()}.txt > ${segments.getSimpleName()}_scarHRD_input.txt 
    """
    """
    #!/usr/bin/env Rscript
    library(scarHRD)
    scarHRD_input <- "${segments.getSimpleName()}_scarHRD_input.txt"
    scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=TRUE)
    """
}


process determine_sex {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    errorStrategy 'retry'

    cpus 1
    memory "2 GB"

    input:
        tuple path(normal_bam), path(normal_bam_bai), path(tumor_bam), path(tumor_bam_bai)

    output:
        tuple val("${normal_bam.getSimpleName()}"), path(result)

    script:

    """
    bam=${normal_bam}

    x_mapped=\$(samtools idxstats \$bam | grep -wE "chrX|X" | cut -f 3)
    x_sequence_length=\$(samtools idxstats \$bam | grep -wE "chrX|X" | cut -f 2)

    y_mapped=\$(samtools idxstats \$bam | grep -wE "chrY|Y" | cut -f 3)
    y_sequence_length=\$(samtools idxstats \$bam | grep -wE "chrY|Y" | cut -f 2)

    xcov=\$(echo "scale=4; \$x_mapped/\$x_sequence_length" | bc)
    ycov=\$(echo "scale=4; \$y_mapped/\$y_sequence_length" | bc)

    if [[ \$ycov == 0 ]]; then
        echo 'female' > result
    else
      ratio=\$(echo "scale=4; \$xcov/\$ycov" | bc)
      echo "X:Y ratio: \$ratio"

      if [[ \$ratio -ge 2 ]]; then
          echo 'female' > result
      else
          echo 'male' > result
      fi
    fi

    exit 0
    """
}

process generate_sequenza_wig {
  conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6"
  errorStrategy 'retry'

  input:
      path(refgenome)

  output:
     path("${refgenome}.gz"), emit: fasta_gz 
     path("${refgenome}_gc50.wig.gz"), emit: wig

  script:
  // 200 for wgs
  // 50 for wes
  def windowsize = 50
  """
  bgzip -c ${refgenome} > ${refgenome}.gz
  sequenza-utils gc_wiggle --fasta ${refgenome}.gz -w $windowsize -o "${refgenome}_gc50.wig.gz"
  """
}

workflow sequenza {

    take:
       matched_bams
       refgenome
       scar_hrd_header

   main:
       gc_wig_file = generate_sequenza_wig(refgenome)
       out = prepare_bam_for_sequenza(matched_bams,gc_wig_file.fasta_gz,gc_wig_file.wig)
       //out.view()
       merge_seqz_files(out.seqzfile)
       formatting_merged_seqz_file(merge_seqz_files.out.header_seqzfile,merge_seqz_files.out.merged_seqzfile)
       seqz_binning(formatting_merged_seqz_file.out)
       sex = determine_sex(matched_bams)
       sex.view()
       smallbin_sex = seqz_binning.out.smallbinfile.join(sex)
       smallbin_sex.view()

       seq_out = r_sequenza(smallbin_sex)
       //hrd_header = Channel.value(scar_hrd_header) 
       preprare_scar_hrd_input(seq_out.segments,seq_out.cellularity_ploidy,scar_hrd_header)
       scar_hrd(preprare_scar_hrd_input.out.scarhrd_input)
       //scar_out = scar_hrd(seq_out.segments,seq_out.cellularity_ploidy)

   //emit:
     //  sequenza_results = seq_out.sequenza_results
      // scar_hrd_results = scar_out.scar_hrd_results
}

workflow sequenza2 {

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

	csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, row.female, file(row.normal_bam), file(row.tumor_bam)) }
	ref_fasta_gz_ch = Channel.value(params.ref_fasta_gz)
	gc_wig_file_ch = Channel.value(params.gc_wig_file)
	outdir = params.outdir
	scar_hrd_header_ch = Channel.value(params.scar_hrd_header)

    sequenza2(csv_ch,ref_fasta_gz_ch,gc_wig_file_ch,scar_hrd_header_ch)
}
