nextflow.enable.dsl=2

process extract_id_bam_filepath {
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'
   
    //conda "conda-forge::pandas=2.2.3"
    conda "conda-forge::python=3.9.15 conda-forge::pandas=2.0.3"
    container "biocontainers/mulled-v2-0594c09780adaaa41fe60b1869ba41c8905a0c98:24a8102d6795963b77f04bb83cc82c081e4a2adc-0"
    
    publishDir "${params.outdir}/", mode: "copy", pattern: "cio_samplesheet.csv"

    input:
      path(target_dir)
      
    output:
      path("bam_path.csv")     , emit: id_bam
      path("purple_path.csv")  , emit: id_purple
      path("cio_samplesheet.csv")  , emit: id_snvs
      
    script:
    """
    id_bam_path.py \\
        --target_dir ${target_dir} \\
        --id_bam_path bam_path.csv \\
        --id_purple_path purple_path.csv \\
        --id_snvs_path cio_samplesheet.csv
    """  
}

process bam_sort_index {

    tag "${bam.getSimpleName()}"
    cpus 8
    memory = { Math.max(16, (task.attempt * bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::samtools=1.22.1"

    input:
        tuple val(sample_ID), path(bam), path(bai)

    output:
        tuple val(sample_ID), path("sorted_${bam.getSimpleName()}.bam"), path("sorted_${bam.getSimpleName()}.bam.bai"), emit: si_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o sorted_${bam.getSimpleName()}.bam
    samtools index -@ ${task.cpus} sorted_${bam.getSimpleName()}.bam -o sorted_${bam.getSimpleName()}.bam.bai
    """
}
    
process bam_qc {
    
    tag "${sorted_bam.getSimpleName()}"
    cpus { sorted_bam.size() > 35.GB ? 8 : 16 }
    memory { sorted_bam.size() > 35.GB ? '128 GB' : '32 GB' }
    //memory = { Math.max(16, (task.attempt * sorted_bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true
    cache 'lenient'

    conda "bioconda::samtools=1.22.1 bioconda::sambamba bioconda::mosdepth"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/coverage", mode: "copy"

    input:
        tuple val(sample_ID), path(sorted_bam), path(sorted_bai)
        val(wes_bedfile)
        //val(tmpdir)
        val(library_type)
        val(grz_bedfile)

    output:
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}.samtools.depth"), path("${sorted_bam.getSimpleName()}.mosdepth.region"), path("${sorted_bam.getSimpleName()}.mosdepth.summary"), path("${sorted_bam.getSimpleName()}.bamqc")
        tuple val(sample_ID), path("filtered_${sorted_bam.getSimpleName()}_markdup.bam"), emit: filtered_bam

    script:
    def sample_name = "${sorted_bam.getSimpleName()}"
    //def  normal_pattern = "NB"
    //def cmd1 = (sample_name =~ /[normal_pattern]/) ? awk1 : awk2
    //def cmd2 = (sample_name =~ /[normal_pattern]/) ? awk3 : awk4
    """
    # markdup with sambamba
    mkdir tmp
    sambamba markdup -t ${task.cpus} --tmpdir tmp ${sorted_bam} ${sorted_bam.getSimpleName()}_markdup.bam
    # bam QC with samtools
    if [[ "${sample_name}" =~ [NB] && "${library_type}" == "wes" ]]; then
        cov_value=30
    elif [[ "${sample_name}" =~ [T] && "${library_type}" == "wes" ]]; then
        cov_value=100
    elif [[ "${sample_name}" =~ [NB] && "${library_type}" == "wgs" ]]; then
        cov_value=20
    elif [[ "${sample_name}" =~ [T] && "${library_type}" == "wgs" ]]; then
        cov_value=30
    else
        echo "No matching filename or wrong library_type"
        exit 1
    fi
    if [[ "${library_type}" == "wes" ]]; then
         bedfile=${wes_bedfile}
    elif [[ "${library_type}" == "wgs" ]]; then
         bedfile=${grz_bedfile}
    else
        echo "Wrong library_type"
        exit 1
    fi
   
    # bam QC with samtools depth without duplicates  
    samtools index -@ ${task.cpus} ${sorted_bam.getSimpleName()}_markdup.bam -o ${sorted_bam.getSimpleName()}_markdup.bam.bai
    samtools view -b -F 1024 -@ ${task.cpus} ${sorted_bam.getSimpleName()}_markdup.bam > filtered_${sorted_bam.getSimpleName()}_markdup.bam
    if [[ "${library_type}" == "wes" ]]; then
        samtools depth -b \$bedfile -aa -@ ${task.cpus} -s filtered_${sorted_bam.getSimpleName()}_markdup.bam > filtered_${sorted_bam.getSimpleName()}_markdup.depth
    elif [[ "${library_type}" == "wgs" ]]; then
        samtools depth -aa -@ ${task.cpus} -s filtered_${sorted_bam.getSimpleName()}_markdup.bam > filtered_${sorted_bam.getSimpleName()}_markdup.depth
    else
        echo "Wrong library_type"
        exit 1
    fi
    cat  filtered_${sorted_bam.getSimpleName()}_markdup.depth | awk -v OFS=',' -v threshold=\$cov_value '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > filtered_${sorted_bam.getSimpleName()}_markdup.cov
    # prepare output samtools depth
    echo "file,mean_cov,targets_above_mincov" > header.csv
    sam_results=\$(cat filtered_${sorted_bam.getSimpleName()}_markdup.cov)
    sam_file=\$(echo "samtools_depth")
    sam_file_results=\$(echo \$sam_file,\$sam_results)
    echo \$sam_file_results > sam_file_results.csv
    cat header.csv sam_file_results.csv > ${sorted_bam.getSimpleName()}.samtools.depth

    # bam QC with samtools depth including duplicated (old)
    if [[ "${library_type}" == "wes" ]]; then
        samtools depth -b \$bedfile -@ ${task.cpus} ${sorted_bam.getSimpleName()}_markdup.bam > with_duplicates_${sorted_bam.getSimpleName()}_markdup.depth
    elif [[ "${library_type}" == "wgs" ]]; then
        samtools depth -@ ${task.cpus} ${sorted_bam.getSimpleName()}_markdup.bam > with_duplicates_${sorted_bam.getSimpleName()}_markdup.depth
    else
        echo "Wrong library_type"
        exit 1
    fi
    cat with_duplicates_${sorted_bam.getSimpleName()}_markdup.depth | awk -v OFS=',' -v threshold=\$cov_value '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > with_duplicates_${sorted_bam.getSimpleName()}_markdup.cov
    # prepare output samtools depth with duplicates
    sam_results_dup=\$(cat with_duplicates_${sorted_bam.getSimpleName()}_markdup.cov)
    sam_file_dup=\$(echo "samtools_depth_w_dup")
    sam_file_results_dup=\$(echo \$sam_file_dup,\$sam_results_dup)
    echo \$sam_file_results_dup > sam_file_dup_results.csv

    # bam QC with mosdepth (for comparison)
    if [[ "${library_type}" == "wes" ]]; then
        mosdepth --by \$bedfile ${sorted_bam.getSimpleName()}_markdup ${sorted_bam.getSimpleName()}_markdup.bam
        zcat  ${sorted_bam.getSimpleName()}_markdup.regions.bed.gz | awk '{len = \$3 - \$2; sum += len * \$4; total += len} END {print sum/total}' > ${sorted_bam.getSimpleName()}.mosdepth.depth
        zcat  ${sorted_bam.getSimpleName()}_markdup.regions.bed.gz | awk -v mincov=\$cov_value '{len = \$3 - \$2; total += len; if (\$4 >= mincov) covered += len} END {print (covered/total)}' > ${sorted_bam.getSimpleName()}.mosdepth.ontargets
        paste -d',' ${sorted_bam.getSimpleName()}.mosdepth.depth ${sorted_bam.getSimpleName()}.mosdepth.ontargets > ${sorted_bam.getSimpleName()}.mosdepth.region
        cat ${sorted_bam.getSimpleName()}_markdup.mosdepth.summary.txt | awk '\$1=="total_region" {print \$4}' > ${sorted_bam.getSimpleName()}.mosdepth.summary
    elif [[ "${library_type}" == "wgs" ]]; then
        mosdepth ${sorted_bam.getSimpleName()}_markdup ${sorted_bam.getSimpleName()}_markdup.bam
        zcat  ${sorted_bam.getSimpleName()}_markdup.per-base.bed.gz | awk '{len = \$3 - \$2; sum += len * \$4; total += len} END {print sum/total}' > ${sorted_bam.getSimpleName()}.mosdepth.depth
        mosdepth --by \$bedfile ${sorted_bam.getSimpleName()}_markdup_regions ${sorted_bam.getSimpleName()}_markdup.bam
        zcat  ${sorted_bam.getSimpleName()}_markdup_regions.regions.bed.gz | awk -v mincov=\$cov_value '{len = \$3 - \$2; total += len; if (\$4 >= mincov) covered += len} END {print (covered/total)}' > ${sorted_bam.getSimpleName()}.mosdepth.ontargets
        paste -d',' ${sorted_bam.getSimpleName()}.mosdepth.depth ${sorted_bam.getSimpleName()}.mosdepth.ontargets > ${sorted_bam.getSimpleName()}.mosdepth.region
        cat ${sorted_bam.getSimpleName()}_markdup.mosdepth.summary.txt | awk '\$1=="total" {print \$4}' > ${sorted_bam.getSimpleName()}.mosdepth.summary
    else
        echo "Wrong library_type"
        exit 1
    fi
    # prepare output mosdepth region
    mosdepth_region_results=\$(cat ${sorted_bam.getSimpleName()}.mosdepth.region)
    mosdepth_region_file=\$(echo "mosdepth_region")
    mosdepth_region_file_results=\$(echo \$mosdepth_region_file,\$mosdepth_region_results)
    echo \$mosdepth_region_file_results > mosdepth_region_file_results.csv
    # prepare output mosdepth summary
    mosdepth_summary_results=\$(cat ${sorted_bam.getSimpleName()}.mosdepth.summary)
    mosdepth_summary_file=\$(echo "mosdepth_summary")
    mosdepth_summary_file_results=\$(echo \$mosdepth_summary_file,\$mosdepth_summary_results)
    echo \$mosdepth_summary_file_results > mosdepth_summary_file_results.csv
    # join results   
    cat ${sorted_bam.getSimpleName()}.samtools.depth  >> ${sorted_bam.getSimpleName()}.bamqc
    cat sam_file_dup_results.csv >> ${sorted_bam.getSimpleName()}.bamqc
    cat mosdepth_region_file_results.csv >> ${sorted_bam.getSimpleName()}.bamqc
    cat mosdepth_summary_file_results.csv >> ${sorted_bam.getSimpleName()}.bamqc 
    """
}

process index_filtered_bams {

    tag "${sample_ID}"
    cpus 8
    memory = { Math.max(16, (task.attempt * filtered_tumorbam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::samtools=1.22.1"

    input:
        tuple val(sample_ID), path(filtered_normalbam), path(filtered_tumorbam)

    output:
        tuple val(sample_ID), path("${filtered_normalbam}"), path("${filtered_normalbam}.bai"), path("${filtered_tumorbam}"), path("${filtered_tumorbam}.bai"), emit: idx_filtered_bam

    script:
    """
    samtools index -@ ${task.cpus} ${filtered_normalbam} -o ${filtered_normalbam}.bai
    samtools index -@ ${task.cpus} ${filtered_tumorbam} -o ${filtered_tumorbam}.bai
    """ 

}

process msisensor_pro {

    tag "${sample_ID}"
    //cpus 16
    memory = { Math.max(16, (task.attempt * idx_filtered_tumorbam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::msisensor-pro=1.2.0"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/msisensor_pro", mode: "copy"

    input:
        tuple val(sample_ID), path(idx_filtered_normalbam), path(idx_filtered_normalbai), path(idx_filtered_tumorbam), path(idx_filtered_tumorbai) 
        val(ref_genome_fasta)
        //val(ref_genome_fai)
        //val(ref_genome_dict)

    output:
        tuple val(sample_ID), path("${idx_filtered_normalbam.getSimpleName()}_${idx_filtered_tumorbam.getSimpleName()}_msi.csv")
        //path("${ref_genome_fasta}_msi.list")
    script:
    """
    msisensor-pro scan \\
        -d ${ref_genome_fasta} \\
        -o msi.list

    msisensor-pro msi \\
        -d msi.list \\
        -n ${idx_filtered_normalbam} \\
        -t ${idx_filtered_tumorbam} \\
        -o ${idx_filtered_normalbam.getSimpleName()}_${idx_filtered_tumorbam.getSimpleName()}_msi.csv
    """
}

process determine_sex {

    tag "${sample_ID}"
    conda "bioconda::samtools=1.22.1"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    errorStrategy 'retry'

    cpus 1
    memory "2 GB"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/coverage", mode: "copy"

    input:
        tuple val(sample_ID), path(normal_bam), path(normal_bam_bai), path(tumor_bam), path(tumor_bam_bai)

    output:
        tuple val(sample_ID), path("${tumor_bam.getSimpleName()}_result_determine_sex.txt"), emit: gender

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
        echo 'female' > ${tumor_bam.getSimpleName()}_result_determine_sex.txt
    else
      ratio=\$(echo "scale=4; \$xcov/\$ycov" | bc)
      echo "X:Y ratio: \$ratio"

      if [[ \$ratio -ge 2 ]]; then
          echo 'female' > ${tumor_bam.getSimpleName()}_result_determine_sex.txt
      else
          echo 'male' > ${tumor_bam.getSimpleName()}_result_determine_sex.txt
      fi
    fi

    exit 0
    """
}

process extract_purple_gender {

    tag "${sample_ID}"
    errorStrategy 'retry'

    cpus 1
    memory "2 GB"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/coverage", mode: "copy"

    input:
        tuple val(sample_ID), path(purple_qc)

    output:
        tuple val(sample_ID), path("${purple_qc.getSimpleName()}_result_purple_sex.txt"), emit: purple_gender
        
    script:
    """
    ambergender=\$(awk -F'\t' '\$1=="AmberGender" {print \$2}' ${purple_qc})
    cobaltgender=\$(awk -F'\t' '\$1=="CobaltGender" {print \$2}' ${purple_qc})

    if [[ "\$ambergender" == "\$cobaltgender" ]]; then
        if [[ "\$ambergender" == "FEMALE" ]]; then
            echo "female" > ${purple_qc.getSimpleName()}_result_purple_sex.txt           
        elif [[ "\$ambergender" == "MALE" ]]; then
            echo "male" > ${purple_qc.getSimpleName()}_result_purple_sex.txt            
        else
            echo "\$ambergender is neither female nor male!"
            exit 1
        fi
    else
        echo "\$ambergender does not match \$cobaltgender"
        exit 1
    fi
    """
}

// ascat obtained from https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf and modified

process ascat {
    tag "${sample_ID}"
    cpus 8

    conda "bioconda::ascat=3.2.0 bioconda::cancerit-allelecount=4.3.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cf02c7911ee5e974ce7db978810770efbd8d872ff5ab3462d2a11bcf022fab5/data'
        : 'community.wave.seqera.io/library/ascat_cancerit-allelecount:c3e8749fa4af0e99'}"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/ascat_cnv", mode: "copy"

    input:
    tuple val(sample_ID), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor), path(purple_gender)
    path allele_files
    path loci_files
    path bed_file
    path fasta
    path gc_file
    path rt_file

    output:
    tuple val(sample_ID), path("*alleleFrequencies_chr*.txt"), emit: allelefreqs
    tuple val(sample_ID), path("*BAF.txt"),                    emit: bafs
    tuple val(sample_ID), path("*cnvs.txt"),                   emit: cnvs
    tuple val(sample_ID), path("*LogR.txt"),                   emit: logrs
    tuple val(sample_ID), path("*metrics.txt"),                emit: metrics
    tuple val(sample_ID), path("*png"),                        emit: png
    tuple val(sample_ID), path("*purityploidy.txt"),           emit: purityploidy
    tuple val(sample_ID), path("*segments.txt"),               emit: segments
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${input_tumor.getSimpleName()}"
    //def gender_value = file("${gender}").readLines().first()
    def gender_file = "${purple_gender}"
    //def gender_value = file(gender_file).readLines().first()
    def genomeVersion = "hg38"
    def purity = "NULL"
    def ploidy = "NULL"
    def gc_input = gc_file            ? "${gc_file}"               : "NULL"
    def rt_input = rt_file            ? "${rt_file}"               : "NULL"
    def minCounts = 20
    def bed_file_arg  = bed_file      ? "BED_file = '${bed_file}',": ""
    def chrom_names_arg = ""
    def min_base_qual = 20
    def min_map_qual = 35
    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    if(dir.exists("${allele_files}")) {
        # expected production use of a directory
        allele_path   = normalizePath("${allele_files}")
        allele_prefix = paste0(allele_path, "/", "${allele_files}", "_chr")
    } else if(file.exists("${allele_files}")) {
        # expected testing use of a single file
        allele_path   = basename(normalizePath("${allele_files}"))
        allele_prefix = sub('_chr[0-9]+\\\\.txt\$', "_chr", allele_path)
    } else {
        stop("The specified allele files do not exist.")
    }

    if(length(Sys.glob(paste0(allele_prefix,"*")) ) == 0) {
        stop(paste("No allele files found matching", allele_prefix))
    }

    if(dir.exists("${loci_files}")) {
        # expected production use of a directory
        loci_path   = normalizePath("${loci_files}")
        loci_prefix = paste0(loci_path, "/", "${loci_files}", "_chr")
    } else if(file.exists("${loci_files}")) {
        # expected testing use of a single file
        loci_path   = basename(normalizePath("${loci_files}"))
        loci_prefix = sub('_chr[0-9]+\\\\.txt\$', "_chr", loci_path)
    } else {
        stop("The specified loci files do not exist.")
    }

    if(length(Sys.glob(paste0(loci_prefix,"*")) ) == 0) {
        stop(paste("No loci files found matching", loci_prefix))
    }

   # Extract gender
   sex <- read.table(file="${gender_file}", header=FALSE)
   if (length(sex) == 1) {
      if (sex == "female") {
        gender_value <- "XX"
      } else if ( sex == "male") {
        gender_value <- "XY"
      } else {
        stop("File contains unexpected content.")
      }
    } else {
      stop("File does not contain exactly one line.")
    }

   # Prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "${input_tumor}",
        normalseqfile = "${input_normal}",
        tumourname = paste0("${prefix}", ".tumour"),
        normalname = paste0("${prefix}", ".normal"),
        allelecounter_exe = "alleleCounter",
        alleles.prefix = allele_prefix,
        loci.prefix = loci_prefix,
        gender = gender_value,
        genomeVersion = "${genomeVersion}",
        nthreads = ${task.cpus},
        minCounts = ${minCounts},
        min_base_qual = ${min_base_qual},
        min_map_qual = ${min_map_qual},
        ${bed_file_arg}
        seed = 42
    )

    # Load the data
    ascat.bc = ascat.loadData(
        Tumor_LogR_file = paste0("${prefix}", ".tumour_tumourLogR.txt"),
        Tumor_BAF_file = paste0("${prefix}", ".tumour_tumourBAF.txt"),
        Germline_LogR_file = paste0("${prefix}", ".tumour_normalLogR.txt"),
        Germline_BAF_file = paste0("${prefix}", ".tumour_normalBAF.txt"),
        genomeVersion = "${genomeVersion}",
        gender = gender_value
    )

    # Plot the raw data
    ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".before_correction."))

    # Optional LogRCorrection
    if("${gc_input}" != "NULL") {

        if(dir.exists("${gc_input}")) {
            # sarek production use of an unzipped folder containing one file
            gc_input = list.files("${gc_input}", recursive = TRUE, full.names = TRUE)
            if(length(gc_input) != 1 | !file.exists(gc_input)) {
                stop("A single gc_input should be provided!")
            }
        } else if(file.exists("${gc_input}")) {
            gc_input = normalizePath("${gc_input}")
        } else {
            stop("gc_input must be a file or folder containing one file")
        }

        if("${rt_input}" != "NULL"){

            if(dir.exists("${rt_input}")) {
                # sarek production use of an unzipped folder containing one file
                rt_input = list.files("${rt_input}", recursive = TRUE, full.names = TRUE)
                if(length(rt_input) != 1 | !file.exists(rt_input)) {
                    stop("A single rt_input should be provided!")
                }
            } else if(file.exists("${rt_input}")) {
                rt_input = normalizePath("${rt_input}")
            } else {
                stop("rt_input must be a file or folder containing one file")
            }

            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = rt_input)
            # Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".after_correction_gc_rt."))
        }
        else {
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = ${rt_input})
            # Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("${prefix}", ".after_correction_gc."))
        }
    }

    # Segment the data
    ascat.bc = ascat.aspcf(ascat.bc, seed=42)

    # Plot the segmented data
    ascat.plotSegmentedData(ascat.bc)

    # Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination,
    # and discrete copy numbers
    # If psi and rho are manually set:
    if (!is.null(${purity}) && !is.null(${ploidy})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=${purity}, psi_manual=${ploidy})
    } else if(!is.null(${purity}) && is.null(${ploidy})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=${purity})
    } else if(!is.null(${ploidy}) && is.null(${purity})){
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=${ploidy})
    } else {
        ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
    }

    # Extract metrics from ASCAT profiles
    QC = ascat.metrics(ascat.bc,ascat.output)

    # Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("${prefix}", ".segments.txt"), sep="\t", quote=F, row.names=F)

    # Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("${prefix}",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    # Write out purity and ploidy info
    summary <- tryCatch({
            matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
                # error handler picks up where error was generated
                print(paste("Could not find optimal solution:  ",err))
                return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
        }
    )
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("${prefix}",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    write.table(QC, file=paste0("${prefix}", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # Version export
    f <- file("versions.yml","w")
    alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
    ascat_version = as.character(packageVersion('ASCAT'))
    writeLines(paste0('"', "${task.process}", '"', ":"), f)
    writeLines(paste("    ascat:", ascat_version), f)
    writeLines(paste("    alleleCounter:", alleleCounter_version), f)
    close(f)
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_ID}"
    """
    touch ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    touch ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    touch ${prefix}.before_correction.test.tumour.germline.png
    touch ${prefix}.before_correction.test.tumour.tumour.png
    touch ${prefix}.cnvs.txt
    touch ${prefix}.metrics.txt
    touch ${prefix}.normal_alleleFrequencies_chr21.txt
    touch ${prefix}.normal_alleleFrequencies_chr22.txt
    touch ${prefix}.purityploidy.txt
    touch ${prefix}.segments.txt
    touch ${prefix}.tumour.ASPCF.png
    touch ${prefix}.tumour.sunrise.png
    touch ${prefix}.tumour_alleleFrequencies_chr21.txt
    touch ${prefix}.tumour_alleleFrequencies_chr22.txt
    touch ${prefix}.tumour_normalBAF.txt
    touch ${prefix}.tumour_normalLogR.txt
    touch ${prefix}.tumour_tumourBAF.txt
    touch ${prefix}.tumour_tumourLogR.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-ascat: \$(Rscript -e "library(ASCAT); cat(as.character(packageVersion('ASCAT')))")
        alleleCounter: \$(alleleCounter --version)
    END_VERSIONS
    """
}

process scarhrd {

    tag "${sample_ID}"

    conda "bioconda::r-sequenza=3.0.0=r41h3342da4_4 conda-forge::r-devtools=2.4.5=r41hc72bb7e_1"

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/scarhrd", mode: "copy"

    input:
      tuple val(sample_ID), path(segments), path(purityploidy)

    output:
      path("${segments.getSimpleName()}.tumour_HRDresults.txt"), emit: scarhrd_results

    script:
    """
    #!/usr/bin/env Rscript
    install.packages('data.table', repos = 'https://cloud.r-project.org')
    library(devtools)
    install_github('sztup/scarHRD',build_vignettes = TRUE)
    install_github('aroneklund/copynumber')

    library(data.table)
    library(scarHRD)
    # prepare scarHRD input
    # get segments and ploidy values
    dt1 <- fread("${segments}", header = TRUE, sep = "\t")
    dt2 <- fread("${purityploidy}", header = TRUE, sep = "\t")
    # extract ploidy
    ploidyValue <- sprintf("%.1f", dt2[[2]][1])
    # add ploidy to segments file
    dt1[, ploidyValue := ploidyValue]
    # calulcate total_cn and add to segments file
    dt1[, total_cn := nMajor + nMinor]
    # reorder segments file
    dt3 <- dt1[, .(sample, chr, startpos, endpos, total_cn, nMajor, nMinor, ploidyValue)]
    # set new column names for scarHRD
    setnames(dt3, c("SampleID", "Chromosome", "Start_position", "End_position", "total_cn", "A_cn", "B_cn", "ploidy"))
    # write scarHRD input file
    fwrite(dt3, file = "scarhrd_input.csv", sep = "\t", quote = FALSE, col.names = TRUE)
    # run scarHRD
    scarHRD_input <- "scarhrd_input.csv"
    scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=FALSE)
    """
}

process format_cnv_table {
    conda "conda-forge::r-base"

    input:
    tuple val(sample_ID), path(segments), path(plurityploidy)

    output:
    tuple val(sample_ID), path("formatted_${segments}")

    script:
    """
    Rscript - <<'RSCRIPT'
    seg <- read.delim("${segments}", check.names = FALSE)

    out <- data.frame(
      SampleID       = seg\$sample,
      Chromosome     = sub("^chr", "", seg\$chr),
      Start_position = seg\$start,
      End_position   = seg\$end,
      total_cn       = seg\$nMajor + seg\$nMinor,
      A_cn           = seg\$nMajor,
      B_cn           = seg\$nMinor
    )

    write.table(
      out,
      "formatted_${segments}",
      sep = "\\t",
      quote = FALSE,
      row.names = FALSE
    )
    RSCRIPT
    """
}

process mini_scarhrd {

    tag "${sample_ID}"
    conda "r-scarhrdcore=0.2.0=r46_0"

    //conda "bioconda::r-sequenza=3.0.0=r41h3342da4_4 conda-forge::r-devtools=2.4.5=r41hc72bb7e_1"
    // container "scarhrdcore"
    containerOptions '--entrypoint /bin/bash'
    conda = null

    publishDir "${params.outdir}/${params.library_type}_${sample_ID}/scarhrd", mode: "copy"

    input:
      // tuple val(sample_ID), path(segments)
      tuple val(sample_ID), path(segments), path(plurityploidy)

    output:
      path("${segments.getSimpleName()}.tumour_HRDresults.txt"), emit: scarhrd_results

    script:
    """
    hrd-score --input $segments --format ascat --assembly grch38 --output "${segments.getSimpleName()}.tumour_HRDresults.txt"
    """
}

workflow postprocessing_wes_custom_oa {

    take:
       target_dir_ch
       wes_bedfile_ch
       //tmpdir_ch
       ref_genome_fasta_ch
       ref_genome_fai_ch
       ref_genome_dict_ch
       library_type_ch
       outdir_ch
       grz_bedfile_ch 
       wes_allele_files_ch
       wes_loci_files_ch
       wes_gc_file
       wes_rt_file

    main:

       extract_id_bam_filepath(target_dir_ch)

       // channel for process_bamfile
       id_bam_qc_ch = extract_id_bam_filepath.out.id_bam | splitCsv(header: true) | map { row -> 
                       [[row.patient_id,
                           file(row.normal_bam),
                           file(row.normal_bai),
                           ],
                       [row.patient_id,
                           file(row.tumor_bam),
                           file(row.tumor_bai),
                           ]
                       ]
                       } | flatMap{ it -> [it[0], it[1]] }

       id_bam_qc_ch.view()

       // channel purple
       id_purple_qc_ch = extract_id_bam_filepath.out.id_purple | splitCsv(header: true) | map { row -> [row.patient_id,file(row.purple_qc)]}

       // channel cio snvs
       id_snvs_ch = extract_id_bam_filepath.out.id_snvs | splitCsv(header: true) | map { row -> [row.sample,file(row.vcf)]}

       // BAM QC
       sorted_and_indexed_bam = bam_sort_index(id_bam_qc_ch)
       bam_qc(sorted_and_indexed_bam.si_bam,wes_bedfile_ch,library_type_ch,grz_bedfile_ch)
       //bam_qc.out.filtered_bam.view() 
       // join and sort tumor normal bam

       // Prepare for msisensor pro
       joined_bams = bam_qc.out.filtered_bam.groupTuple() 
       //joined_bams.view()
       sorted_joined_bams = joined_bams.map { id, files ->
                                           def order = ['N': 0, 'B': 1, 'T': 2]
                                           //tuple(id, *files.sort { order[(it.name =~ /[NBT]/)[0]] })
                                           tuple(id, *files.sort { f ->
                                               //def letter = (f.name =~ /[NBT]/)[0]
                                               def letter = f.name.find(/[NBT]/)
                                               order[letter]
                                            })
                                        }
       //sorted_joined_bams.view()
       idx_filtered_bam = index_filtered_bams(sorted_joined_bams)
       idx_filtered_bam.view()

       // Msisensor pro
       // msisensor_pro(idx_filtered_bam,ref_genome_fasta_ch)

       // Determine sex has to be fixed !
       // gender = determine_sex(idx_filtered_bam)
       //gender.view()

       // Extract gender from purple_qc
       gender_purple = extract_purple_gender(id_purple_qc_ch)

       // Ascat and scarHRD
       idx_filtered_bam_gender = idx_filtered_bam.join(gender_purple, by:0)
       idx_filtered_bam_gender.view()
       ascat(idx_filtered_bam_gender,wes_allele_files_ch,wes_loci_files_ch,wes_bedfile_ch,ref_genome_fasta_ch,wes_gc_file,wes_rt_file)
       segments_purityploidy = ascat.out.segments.join(ascat.out.purityploidy, by:0)
       //segments_purityploidy.view()

       /* 
       old scarhrd
       scarhrd(segments_purityploidy)                                                  
       */

       // mini_scarhrd(format_cnv_table(segments_purityploidy))
       mini_scarhrd(segments_purityploidy)
}

workflow {

    // set channels
    target_dir_ch = Channel.fromPath(params.target_dir, type: "dir")
    wes_bedfile_ch = Channel.value(params.wes_bedfile)
    //tmpdir_ch = Channel.value(params.tmpdir)
    ref_genome_fasta_ch = Channel.value(params.ref_genome_fasta)
    ref_genome_fai_ch = Channel.value(params.ref_genome_fai)
    ref_genome_dict_ch = Channel.value(params.ref_genome_dict)
    library_type_ch = Channel.value(params.library_type)
    outdir_ch = Channel.value(params.outdir)
    grz_bedfile_ch = Channel.value(params.grz_bedfile)
    wes_allele_files_ch = Channel.value(params.wes_allele_files)
    wes_loci_files_ch = Channel.value(params.wes_loci_files)
    wes_gc_file = Channel.value(params.wes_gc_file)
    wes_rt_file = Channel.value(params.wes_rt_file)

    postprocessing_wes_custom_oa(target_dir_ch,
                                 wes_bedfile_ch,
                                 ref_genome_fasta_ch,
                                 ref_genome_fai_ch,
                                 ref_genome_dict_ch,
                                 library_type_ch,
                                 outdir_ch,
                                 grz_bedfile_ch,
                                 wes_allele_files_ch,
                                 wes_loci_files_ch,
                                 wes_gc_file,
                                 wes_rt_file)
}
