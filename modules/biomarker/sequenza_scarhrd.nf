nextflow.enable.dsl=2

process bam_proc {

    tag "${sample_id}"
    debug true
            
    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7' 

    cpus 24
    errorStrategy 'retry'
    memory '250 GB'

    input:
        tuple path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
        path(ref_fasta)
        path(gc_wig_file)
        
    output:
        tuple val(sample_id), path("${tumor_bam.getSimpleName()}*.seqz.gz"), emit: seqzfiles
        tuple val(sample_id), path("${tumor_bam.getSimpleName()}*.seqz.gz.tbi"), emit: tbi_seqzfiles
            
    script:
    def args = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    sample_id = "${tumor_bam.getSimpleName()}"
    """
    samtools faidx ${ref_fasta}
   
    sequenza-utils bam2seqz \\
                   --normal ${normal_bam} \\
                   --tumor ${tumor_bam} \\
                   --fasta ${ref_fasta} \\
                   -C $args \\
                   --parallel ${task.cpus} \\
                   -gc ${gc_wig_file} \\
                   --output ${tumor_bam.getSimpleName()}.seqz.gz
    """
}

process seqz_proc {

    tag "${sample_id}"
    debug true

    conda "bioconda::sequenza-utils=3.0.0=py39he10ea66_6"
    //container 'biocontainers/sequenza-utils:3.0.0--py312h719dbc0_7'

    cpus 24
    errorStrategy 'retry'

    input:
        tuple val(sample_id), path(seqzfiles)

    output:
        tuple val(sample_id), path("${sample_id}_bin.seqz.gz"), emit: smallbinfile

    script:
    def header_seqz="chromosome      position        base.ref        depth.normal    depth.tumor     depth.ratio     Af      Bf      zygosity.normal GC.percent      good.reads      AB.normal   AB.tumor tumor.strand"
    """
    for file in `ls -v ${seqzfiles}` 
    do 
        zcat \${file} | awk 'NR>1 {print}' >> merged_${sample_id}.seqz
    done 
     
    echo $header_seqz > header_${sample_id}.txt
    cat header_${sample_id}.txt merged_${sample_id}.seqz | bgzip -@ ${task.cpus} > all_${sample_id}.seqz.gz  
    tabix -f -s 1 -b 2 -e 2 -S 1 all_${sample_id}.seqz.gz

    sequenza-utils seqz_binning \\
              --seqz all_${sample_id}.seqz.gz \\
              --window 50 \\
              -o ${sample_id}_bin.seqz.gz
    """
}

process r_sequenza {

    tag "${sample_id}"
    debug true

    conda "/home/pbasitta/miniconda3/envs/r_seq"

    errorStrategy 'retry'

    input:
       tuple val(sample_id), path(smallbinfile), val(female)

    output:
       tuple val(sample_id), path("${sample_id}_segments.txt"), path("${sample_id}_cellularity_ploidy.txt"), emit: segments_ploidy
       path("${sample_id}_sequenza_results.tar.gz"), emit: sequenza_results
       //tuple val(sample_id), path("${sample_id}_HRDresults.txt"), emit: scar_hrd_results

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

    dir.create("${sample_id}_sequenza_results")
    pdf_lst <- list.files(".", "*.pdf\$")
    txt_lst <- list.files(".","*.txt\$")
    RData_lst <- list.files(".","*.RData\$")
    file.copy(pdf_lst,"${sample_id}_sequenza_results")
    file.copy(txt_lst,"${sample_id}_sequenza_results")
    file.copy(RData_lst,"${sample_id}_sequenza_results")
    tar("${sample_id}_sequenza_results.tar.gz", "./${sample_id}_sequenza_results", compression = "gzip", tar = "tar")
    """
}

process scar_hrd {

    tag "${sample_id}"
    debug true

    conda "/home/pbasitta/miniconda3/envs/r_seq"

    input:
      tuple val(sample_id), path(segments), path(cellularity_ploidy)
      val(scar_hrd_header)

    output:
       tuple val(sample_id), path("${sample_id}_HRDresults.txt"), emit: scar_hrd_results

    script:
    """
    awk 'BEGIN {FS=OFS="\t"} {print \$1,\$2,\$3,\$10,\$11,\$12}' ${segments} > tmp1_${sample_id}.txt
    sed -i "s/\$/\t${sample_id}/" tmp1_${sample_id}.txt
    ploidyValue=`awk 'NR>1{print \$2}' ${cellularity_ploidy}`
    sed -i "s/\$/\t\$ploidyValue/" tmp1_${sample_id}.txt
    awk 'BEGIN {FS=OFS="\t"} {print \$7,\$1,\$2,\$3,\$4,\$5,\$6,\$8}' tmp1_${sample_id}.txt > tmp2_${sample_id}.txt
    awk 'NR>1 {print}' tmp2_${sample_id}.txt > tmp3_${sample_id}.txt
    cat ${scar_hrd_header} tmp3_${sample_id}.txt > ${sample_id}_scarHRD_input.txt

    scar_hrd.R "${sample_id}_scarHRD_input.txt"
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
        tuple val("${tumor_bam.getSimpleName()}"), path(result)

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
     //path("${refgenome}.gz"), emit: fasta_gz 
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
       bam_proc(matched_bams, refgenome, gc_wig_file.wig)
       out = seqz_proc(bam_proc.out.seqzfiles)
       sex = determine_sex(matched_bams)
       smallbin_sex = out.smallbinfile.join(sex)
       seq_out = r_sequenza(smallbin_sex)
       scar_out = scar_hrd(seq_out.segments_ploidy,scar_hrd_header)
       
   emit:
       sequenza_results = seq_out.sequenza_results
       scar_hrd_results = scar_out.scar_hrd_results
}
