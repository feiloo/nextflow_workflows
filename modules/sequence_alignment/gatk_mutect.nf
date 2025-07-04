// implement best-practice somatic variant calling
// reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2

include { index_bam; index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

process gatk_mutect {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory {Math.min(196, 196+(64 * (task.attempt-1))).GB}

    input:
        tuple path(normal_bam), path(normal_bam_index), path(tumor_bam), path(tumor_bam_index)
	path(intervals)
	path(panel_of_normals)
	path(panel_of_normals_index)
	path(germline_resource)
	path(germline_resource_index)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      tuple path("${tumor_bam.getSimpleName()}_unfiltered.vcf"), path("${tumor_bam.getSimpleName()}_unfiltered.vcf.stats"), emit: vcf
      path("${tumor_bam.getSimpleName()}_f1r2_data.shards.tar"), emit: f1r2_data

    script:
    def normal_samplename = "${normal_bam.getSimpleName()}"
    """

mkdir -p beds

script="
from itertools import groupby

with open('${intervals}') as f:
    txt = f.readlines()


for k,g in groupby(txt, lambda x: x.split('\t')[0]):
    gr = list(g)
    with open(f'beds/{k}.bed', 'w') as f2:
        f2.writelines(gr)
"

python3 -c "\$script"

    export JAVA_HEAP=8
    export HMM_THREADS=12

    Mutect_chrom () {
    INTERVALS=\$1

    mkdir -p tmp/\$2
    gatk Mutect2 \\
	--java-options "-Djava.io.tmpdir=tmp/\$2 -Xms\${JAVA_HEAP}G -Xmx\${JAVA_HEAP}G" \\
	--native-pair-hmm-threads \${HMM_THREADS} \\
	--tmp-dir tmp/\$2 \\
	--input ${tumor_bam} \\
	--input ${normal_bam} \\
	-normal ${normal_bam.getSimpleName()} \\
	-L \$INTERVALS \\
	-pon ${panel_of_normals} \\
	-germline-resource ${germline_resource} \\
        --reference ${refgenome} \\
	--f1r2-tar-gz "${tumor_bam.getSimpleName()}_f1r2_data_\$2.tar.gz" \\
	--output "${tumor_bam.getSimpleName()}_unfiltered_\$2.vcf"

    }
    export -f Mutect_chrom

    # calculate max jobs for the given resources
    MEMJOBS=\$((${task.memory.getGiga()} / \$JAVA_HEAP))
    CPUJOBS=${task.cpus}
    JOBS=\$((MEMJOBS<CPUJOBS ? MEMJOBS : CPUJOBS))

    ls beds/* | parallel --jobs \$JOBS "Mutect_chrom {} {#}" --pipe

    ls *.vcf > outputs.list
    ls *.stats > outputs_stats.list

    gatk MergeVcfs -I outputs.list -O "${tumor_bam.getSimpleName()}_unfiltered.vcf"
    gatk MergeMutectStats --stats outputs_stats.list -O "${tumor_bam.getSimpleName()}_unfiltered.vcf.stats"
    tar -cf "${tumor_bam.getSimpleName()}_f1r2_data.shards.tar" *_f1r2_data_*

    """
}



process gatk_mutect_single {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory {Math.min(56, 36+(10 * (task.attempt-1))).GB}

    input:
        tuple path(normal_bam), path(normal_bam_index), path(tumor_bam), path(tumor_bam_index)
	path(intervals)
	path(panel_of_normals)
	path(panel_of_normals_index)
	path(germline_resource)
	path(germline_resource_index)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      tuple path("${tumor_bam.getSimpleName()}_unfiltered.vcf"), path("${tumor_bam.getSimpleName()}_unfiltered.vcf.stats"), emit: vcf
      path("${tumor_bam.getSimpleName()}_f1r2_data.tar.gz"), emit: f1r2_data

    script:
    def normal_samplename = "${normal_bam.getSimpleName()}"
    """
    mkdir -p tmp
    gatk Mutect2 \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--native-pair-hmm-threads ${task.cpus} \\
	--tmp-dir tmp \\
	--input ${tumor_bam} \\
	--input ${normal_bam} \\
	-normal ${normal_bam.getSimpleName()} \\
	-L ${intervals} \\
	-pon ${panel_of_normals} \\
	-germline-resource ${germline_resource} \\
        --reference ${refgenome} \\
	--f1r2-tar-gz "${tumor_bam.getSimpleName()}_f1r2_data.tar.gz" \\
	--output "${tumor_bam.getSimpleName()}_unfiltered.vcf"
    """
}

process gatk_getsamplename {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    cpus 1

    input:
        path(sample_bam)

    output:
      val(header_samplename)

    script:
    """
    gatk GetSampleName \\
	--INPUT ${sample_bam}
    """
}

process gatk_learn_readorientationmodel {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory {Math.min(56, 40+(10 * (task.attempt-1))).GB}

    input:
        path(f1r2_data)

    output:
      path("${f1r2_data.getSimpleName()}_model.tar.gz"), emit: orientation_model

    script:
    """
    mkdir -p tmp

    tar -xf ${f1r2_data}
    # note, make sure that the input archive itself included in it
    all_f1r2=\$(ls *_f1r2_*.tar.gz | xargs -i printf -- " -I {} ")

    gatk LearnReadOrientationModel \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	\$all_f1r2 \\
	--output "${f1r2_data.getSimpleName()}_model.tar.gz"
    """
}


process gatk_getpileupsummaries {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory {Math.min(56, 40+(10 * (task.attempt-1))).GB}

    input:
        tuple path(sample_bam), path(sample_bam_index)
        path(genomic_intervals)
	// prepared from like gnomAD with pop-freqs in the info field
        path(variant_frequency_vcf)
        path(variant_frequency_vcf_index)

    output:
      path("${sample_bam.getSimpleName()}_pileup.table"), emit: table

    script:
    """

    mkdir -p tmp
    gatk GetPileupSummaries \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	-I ${sample_bam} \\
	-L ${genomic_intervals} \\
	-V ${variant_frequency_vcf} \\
	-O "${sample_bam.getSimpleName()}_pileup.table"
    """
}

process gatk_calculate_contamination {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory '56 GB'

    input:
        tuple path(normal_pileups), path(tumor_pileups)

    output:
      tuple path("${tumor_pileups.getSimpleName()}_contamination.table"), path("${tumor_pileups.getSimpleName()}_segments.tsv"), emit: contamination_table_and_segments

    script:
    """

    mkdir -p tmp
    gatk CalculateContamination \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--input ${tumor_pileups} \\
	-matched ${normal_pileups} \\
	--tumor-segmentation ${tumor_pileups.getSimpleName()}_segments.tsv \\
	--output "${tumor_pileups.getSimpleName()}_contamination.table"
    """
}

process gatk_filter_calls {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    memory '56 GB'

    input:
        tuple path(sample_vcf), path(sample_vcf_stats), path(contamination_table), path(tumor_segments), path(orientation_model)
	path(refgenome)
	path(refgenome_index)
	path(refgenome_dict)

    output:
      path("${outputfile}"), emit: vcf

    script:

    outputfile = "${sample_vcf.getSimpleName().substring(0, sample_vcf.getSimpleName().length() - 11)}.vcf"

    """
    mkdir -p tmp

    gatk FilterMutectCalls \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--variant ${sample_vcf} \\
	--ob-priors ${orientation_model} \\
	--contamination-table ${contamination_table} \\
	--tumor-segmentation ${tumor_segments} \\
	--output "${outputfile}" \\
	-R "${refgenome}"
    """
}

process create_pon_db {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    cpus 1

    input:
    	path(pon_vcfs)
    	path(refgenome)
    	path(intervals)

    output:
        path("pon_genomics_db")

    script:
    """
    mkdir -p tmp
    gatk GenomicsDBImport \\
        -R ${refgenome} \\
        --genomicsdb-workspace-path pon_genomics_db \\
    """
}

def filename_to_dict(filename){
        // split filename to its extensions
        def filename_extensions = filename.tokenize(".")
        // split to its naming-scheme parts
        def parts = filename_extensions[0].tokenize("_")
        return [sampleid:parts[0], sampletype:parts[1], modality:parts[2], readnr:parts[3], suffix:parts[4]]
}


workflow variant_call {
  take:
    sample_bams
    intervals
    panel_of_normals
    germline_resource
    refgenome
    args

  main:

    // best practices from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
    // How to Call somatic mutations using GATK4 Mutect2

    // key is XXX_T_1
    sample_bams_w_key = sample_bams.map{ it -> ["${it.getSimpleName()}", it]}
    sample_bams_idx_w_key = index_bam(sample_bams).map{ it -> ["${it.getSimpleName()}", it] }
    // SAMPLENAME, [XX.bam, XX.bai]
    sample_bams_w_indices = sample_bams_w_key.join(sample_bams_idx_w_key).map{ it -> [filename_to_dict(it[0]).sampleid, [it[1], it[2]]] }

    // warning, whatever the path type is for bam1 and bam2, 
    // the strings and .endsWith() methods do not work as expected
    def pair_bams_fn = { it ->
        def bam1_p = filename_to_dict(it[1][0][0].name)
        def bam1 = it[1][0][0]
        def bam1_idx = it[1][0][1]
        def bam2_p = filename_to_dict(it[1][1][0].name)
        def bam2 = it[1][1][0]
        def bam2_idx = it[1][1][1]
 
        // mutect takes the normal bam first
        // switch if the order is wrong
        if(bam1_p.sampletype == "T" && bam2_p.sampletype == "N"){
            return [bam2, bam2_idx, bam1, bam1_idx]
        } else if(bam1_p.sampletype == "N" && bam2_p.sampletype == "T"){
            return [bam1, bam1_idx, bam2, bam2_idx]
        } else {
          throw new Exception("error ordering bam pair: ${it}")
        }
 
    }

    bam_pairs = sample_bams_w_indices.groupTuple(by: 0, size:2).map{ it -> pair_bams_fn(it) }

    refgenome_index = index_fasta(args.refgenome).fasta_index
    refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict

    indices = Channel.of(germline_resource, panel_of_normals)
    indices_mid = gatk_indexfeaturefile(indices).known_sites_index

    // panel_of_normals_index = gatk_indexfeaturefile(panel_of_normals).known_sites_index
    panel_of_normals_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(panel_of_normals).getSimpleName()}" }
    germline_resource_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(germline_resource).getSimpleName()}" }


    mut = gatk_mutect(bam_pairs, intervals, 
	    panel_of_normals, panel_of_normals_index,
	    germline_resource,  germline_resource_index,  
	    refgenome, refgenome_index, refgenome_dict)

    vcf = mut.vcf.map{it -> it[0]}
    om = gatk_learn_readorientationmodel(mut.f1r2_data).orientation_model

    // for tumors and for normals
    // todo, check that ew can use the same germline resource for pileup here

    sample_bams_w_indices_no_key = sample_bams_w_indices.map{it -> it[1]}
    all_pileups = gatk_getpileupsummaries(sample_bams_w_indices_no_key, intervals, germline_resource, germline_resource_index).table

    def pair_pileups = { it ->
        def pu1_p = filename_to_dict(it[1][0].name)
        def pu1 = it[1][0]
        def pu2_p = filename_to_dict(it[1][1].name)
        def pu2 = it[1][1]

        if(pu1_p.sampletype == "T" && pu2_p.sampletype == "N" ){
            return [pu2, pu1]
        } else if(pu2_p.sampletype == "T" && pu1_p.sampletype == "N" ){
            return [pu1, pu2]
        } else {
          throw new Exception("error ordering pileups pair: ${it}")
        }

    }

    // group pileups by samplename
    matched_pileups = all_pileups.map{it -> ["${filename_to_dict(it.getSimpleName()).sampleid}", it]}.groupTuple(size: 2, sort: true).map{it -> pair_pileups(it)}

    contamination_table_and_segments = gatk_calculate_contamination(matched_pileups).contamination_table_and_segments
    c_w_key = contamination_table_and_segments.map{it -> ["${it[0].getSimpleName().split('_')[0]}", it]}
    vcf_w_key = mut.vcf.map{it -> ["${it[0].getSimpleName().split('_')[0]}", it]}
    om_w_key = om.map{it -> ["${it.getSimpleName().split('_')[0]}", it]}
    vcf_w_filter_data = vcf_w_key.join(c_w_key).join(om_w_key).map{it -> [it[1][0], it[1][1], it[2][0], it[2][1], it[3]]}

    filtered_vcf = gatk_filter_calls(vcf_w_filter_data, args.refgenome, refgenome_index, refgenome_dict).vcf

  emit:
    vcf = filtered_vcf
}

/*
workflow create_pon {
  take:
    args
  main:
    vcf = filtered_vcf
}
*/
