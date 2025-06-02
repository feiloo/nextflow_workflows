nextflow.enable.dsl=2

process prepare_bam_for_sequenza {
            
    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6 bioconda::r-sequenza=3.0.0 conda-forge::r-devtools=2.4.5"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7' 

    //publishDir(path: "${outdir}/${sample_id}/preprocessed_files", mode: "copy")
    cpus 24
    errorStrategy 'ignore'

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

        tuple val(sample_id), path("${sample_id}_bin.seqz.gz"), emit: smallbinfile
            
    script:
    def args = task.ext.args ?: ""     
    def args2 = task.ext.args2 ?: ""
    sample_id = "${normal_bam.getSimpleName()}"
    """
    #!/usr/bin/env Rscript
    library(devtools)
    install_github('aronenklund/copynumber')
    install_github('sztup/scarHRD', build_vignettes = TRUE)
    """


    """
    conda remove r-vroom

    sequenza-utils bam2seqz \\
                   --normal ${normal_bam} \\
                   --tumor ${tumor_bam} \\
                   --fasta ${ref_fasta_gz} \\
                   --parallel ${task.cpus} \\
                   -gc ${gc_wig_file} \\
                   --output ${normal_bam.getSimpleName()}.seqz.gz

    for file in `ls -v *chr[0-9]*.seqz.gz` 
    do 
        zcat \${file} | awk 'NR==1 {print; exit}' > header_${sample_id}.txt
        zcat \${file} | awk 'NR>1 {print}' >> merged_${sample_id}.seqz
    done 

    cat header${sample_id}.txt merged_${sample_id}.seqz | bgzip -@ ${task.cpus} > all_${sample_id}.seqz.gz  
    tabix -f -s 1 -b 2 -e 2 -S 1 all_${sample_id}.seqz.gz

    sequenza-utils seqz_binning \\
              --seqz all_${sample_id}.seqz.gz \\
              --window 200 \\
              -o ${sample_id}_bin.seqz.gz
    """
}


process r_sequenza { 
    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6 bioconda::r-sequenza=3.0.0 conda-forge::r-devtools=2.4.5"
    //conda "bioconda::r-sequenza=3.0.0=r40h3342da4_3"
    
    //publishDir(path: "${outdir}/${sample_id}/sequenza_results", mode: "copy")
    errorStrategy 'ignore'
    
    input:
       tuple val(sample_id), path(smallbinfile), val(female)
  
    output:
       // tuple val(sample_id), path("${sample_id}_alternative_fit.pdf"),        emit: alternative_fit 
       // tuple val(sample_id), path("${sample_id}_alternative_solutions.txt"),  emit: alternative_solutions
       // tuple val(sample_id), path("${sample_id}_chromosome_depths.pdf"),      emit: chromosome_depths
       // tuple val(sample_id), path("${sample_id}_chromosome_view.pdf"),        emit: chromosome_view
       // tuple val(sample_id), path("${sample_id}_CN_bars.pdf"),                emit: CN_bars
       // tuple val(sample_id), path("${sample_id}_confints_CP.txt"),            emit: confints_CP
       // tuple val(sample_id), path("${sample_id}_CP_contours.pdf"),            emit: CP_contours
       // tuple val(sample_id), path("${sample_id}_gc_plots.pdf"),               emit: gc_plots
       // tuple val(sample_id), path("${sample_id}_genome_view.pdf"),            emit: genome_view
       // tuple val(sample_id), path("${sample_id}_model_fit.pdf"),              emit: model_fit
       // tuple val(sample_id), path("${sample_id}_mutations.txt"),              emit: mutations
       // tuple val(sample_id), path("${sample_id}_sequenza_cp_table.RData"),    emit: sequenza_cp_table
       // tuple val(sample_id), path("${sample_id}_sequenza_extract.RData"),     emit: sequenza_extract
       // tuple val(sample_id), path("${sample_id}_sequenza_log.txt"),           emit: sequenza_log

       //tuple val(sample_id), path("${sample_id}_segments.txt"),               emit: segments
       //tuple val(sample_id), path("${sample_id}_cellularity_ploidy.txt"),     emit: cellularity_ploidy

       tuple val("${smallbinfile.getSimpleName()}"), path("${smallbinfile.getSimpleName()}_segments.txt"), path("${smallbinfile.getSimpleName()}_cellularity_ploidy.txt"),     emit: segments_ploidy
       path("${sample_id}_sequenza_results.tar.gz}"), emit: sequenza_results

    script:

    """
    #!/usr/bin/env Rscript
    library(devtools)
    install_github('aronenklund/copynumber')
    install_github('sztup/scarHRD', build_vignettes = TRUE)
    """

    """
    conda remove r-vroom
    """

    """
    #!/usr/bin/env Rscript
    library(sequenza)
    data.file <- "${smallbinfile}"
    seqzdata <- sequenza.extract(data.file, assembly="hg38")

    if (length(file_content) == 1) {
      if (file_content == "female") {
	female <- TRUE
      } else if (file_content == "male") {
	female <- FALSE
      } else {
	stop("File contains unexpected content.")
      }
    } else {
      stop("File does not contain exactly one line.")
    }

    CP <- sequenza.fit(seqzdata, female="${female}")
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${sample_id}", out.dir = ".")
    
    cint <- get.ci(CP)
    cellularity <- cint[["max.cellularity"]]
    ploidy <- cint[["max.ploidy"]]
    data_CP <- data.frame(cellularity = cellularity, ploidy = ploidy)
    write.table(data_CP, file="${sample_id}_cellularity_ploidy.txt", sep="\t", row.names=FALSE)
    """

    """
    mkdir sequenza_results
    cp *.pdf sequenza_results
    cp *.txt sequenza_results
    cp *.RData sequenza_results
    tar czf ${sample_id}_sequenza_results.tar.gz
    """
}     


process scar_hrd {
    
    tag "${sample_id}"
    debug true
    errorStrategy 'ignore'

    publishDir(path: "${outdir}/${sample_id}/scarHRD_results", mode: "copy")

    input:
      tuple val(sample_id), path(segments), path(cellularity_ploidy)


    output:
      path("${sample_id}_HRDresults.txt"), emit: scar_hrd_results

    script:
    def args = task.ext.args ?: ""
    def scar_hrd_header = 'SampleID        Chromosome      Start_position  End_position    total_cn        A_cn    B_cn    ploidy'

    """
    awk 'BEGIN {FS=OFS="\t"} {print \$1,\$2,\$3,\$10,\$11,\$12}' ${segments} > tmp1_${sample_id}.txt
    sed -i "s/\$/\t${sample_id}/" tmp1_${sample_id}.txt
    ploidyValue=`awk 'NR>1{print \$2}' ${cellularity_ploidy}`
    sed -i "s/\$/\t\$ploidyValue/" tmp1_${sample_id}.txt
    awk 'BEGIN {FS=OFS="\t"} {print \$7,\$1,\$2,\$3,\$4,\$5,\$6,\$8}' tmp1_${sample_id}.txt > tmp2_${sample_id}.txt
    awk 'NR>1 {print}' tmp2_${sample_id}.txt > tmp3_${sample_id}.txt

    cat ${scar_hrd_header} tmp3_${sample_id}.txt > ${sample_id}_scarHRD_input.txt 
    """

    """
    #!/usr/bin/env Rscript
    library(scarHRD)
    scarHRD_input <- "${scarhrd_input}"
    scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=$args)
    """
}


process determine_sex {
    conda "bioconda::samtools=1.17"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    errorStrategy 'ignore'

    cpus 1
    memory "2 GB"

    input:
        tuple path(normal_bam), path(normal_bam_bai), path(tumor_bam), path(tumor_bam_bai)

    output:
        tuple val("${normal_bam.getSimpleName()}"), path(result)

    script:

    """
    run_stub=1

    if [[ \$runstub == 1 ]]; then
      echo 'female' > result
      exit 0
    fi

    bam=${normal_bam}

    x_mapped=\$(samtools idxstats \$bam | grep "X" | cut -f 3)
    x_sequence_length=\$(samtools idxstats \$bam | grep "X" | cut -f 2)

    y_mapped=\$(samtools idxstats \$bam | grep "Y" | cut -f 3)
    y_sequence_length=\$(samtools idxstats \$bam | grep "Y" | cut -f 2)

    xcov=\$(echo "scale=4; \$x_mapped/\$x_sequence_length | bc)
    ycov=\$(echo "scale=4; \$y_mapped/\$y_sequence_length | bc)

    if [[ \$ycov == 0 ]]; then
        echo 'female' > result
    else
      rat=\$(echo "scale=4; \${xcov}/\${ycov}" | bc)
      echo "X:Y ratio: \$rat"

      if [[ \$rat >= 4 ]]; then
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
  errorStrategy 'ignore'

  input:
      path(refgenome)

  output:
      path("${refgenome}_gc200.wig.gz"), emit: wig

  script:
  // 200 for wgs
  // 50 for wes
  def windowsize = 200
  """
  sequenza-utils gc_wiggle --fasta ${refgenome} -w 200 -o "${refgenome}_gc200.wig.gz"
  """
}

workflow sequenza {

    take:
       matched_bams
       refgenome

       
   main:
       gc_wig_file = generate_sequenza_wig(refgenome)
       out = prepare_bam_for_sequenza(matched_bams, refgenome, gc_wig_file.wig)
       sex = determine_sex(matched_bams)
       smallbin_sex = out.smallbinfile.join(sex)

       seq_out = r_sequenza(smallbin_sex)
       scar_out = scar_hrd(seq_out.segments_ploidy)

   emit:
       sequenza_results = seq_out.sequenza_results
       scar_hrd_results = scar_out.scar_hrd_results
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
