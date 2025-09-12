// implement best-practice somatic variant calling
// reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2

include { index_bam; index_fasta } from "$NEXTFLOW_MODULES/sequence_alignment/samtools.nf"
include { gatk_indexfeaturefile; gatk_createsequencedictionary } from "$NEXTFLOW_MODULES/sequence_alignment/gatk.nf"

def file_id(f) {
    return f.getSimpleName().split('_')[0..3].join('_')
}
def vcf_name_from_pair(normal_file,tumor_file){
	// assert normal_file.getSimpleName() ==~ "/
	return "${file_id(normal_file)}__${file_id(tumor_file)}"
}

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
      //tuple path("${tumor_bam.getSimpleName()}_unfiltered.vcf"), path("${tumor_bam.getSimpleName()}_unfiltered.vcf.stats"), emit: vcf
      tuple path("${vcf_name_from_pair(normal_bam, tumor_bam)}_unfiltered.vcf"), path("${vcf_name_from_pair(normal_bam,tumor_bam)}_unfiltered.vcf.stats"), emit: vcf
      path("${vcf_name_from_pair(normal_bam, tumor_bam)}_f1r2_data.shards.tar"), emit: f1r2_data

    script:
    def out_prefix = "${vcf_name_from_pair(normal_bam, tumor_bam)}"
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
	--f1r2-tar-gz "${out_prefix}_f1r2_data_\$2.tar.gz" \\
	--output "${out_prefix}_unfiltered_\$2.vcf"

    }
    export -f Mutect_chrom

    # calculate max jobs for the given resources
    MEMJOBS=\$((${task.memory.getGiga()} / \$JAVA_HEAP))
    CPUJOBS=${task.cpus}
    JOBS=\$((MEMJOBS<CPUJOBS ? MEMJOBS : CPUJOBS))

    ls beds/* | parallel --jobs \$JOBS "Mutect_chrom {} {#}" --pipe

    ls *.vcf > outputs.list
    ls *.stats > outputs_stats.list

    gatk MergeVcfs -I outputs.list -O "${out_prefix}_unfiltered.vcf"
    gatk MergeMutectStats --stats outputs_stats.list -O "${out_prefix}_unfiltered.vcf.stats"
    tar -cf "${out_prefix}_f1r2_data.shards.tar" *_f1r2_data_*

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
      tuple path("${vcf_name_from_pair(normal_pileups, tumor_pileups)}_contamination.table"), path("${vcf_name_from_pair(normal_pileups, tumor_pileups)}_segments.tsv"), emit: contamination_table_and_segments

    script:
    def out_prefix = vcf_name_from_pair(normal_pileups, tumor_pileups)
    """
    mkdir -p tmp
    gatk CalculateContamination \\
	--java-options "-Djava.io.tmpdir=tmp -Xms50G -Xmx50G" \\
	--input ${tumor_pileups} \\
	-matched ${normal_pileups} \\
	--tumor-segmentation ${out_prefix}_segments.tsv \\
	--output "${out_prefix}_contamination.table"
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
        return [sampleid:parts[0], sampletype:parts[1], modality:parts[2], material_type:parts[3], readnr:parts[4], suffix:parts[5]]
}


def removeSuffix(String str, String suffix) {
    // this is strict now
    assert str != null
    assert suffix != null
    assert str.endsWith(suffix)
    return str[0..<str.length()-suffix.length()]
}

def create_interval_shards(intervals_file) {
    def chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                      'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                      'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
    
    def chromosomes_data = [:]
    def unknown_chroms = []
    
    intervals_file.eachLine { line ->
        if (!line.startsWith('#') && line.trim()) {
            def parts = line.split('\t')
            def chrom = parts[0]
            if (!chromosomes.contains(chrom)) {
                unknown_chroms.add(chrom)
            }
            if (!chromosomes_data.containsKey(chrom)) {
                chromosomes_data[chrom] = []
            }
            chromosomes_data[chrom].add(line)
        }
    }
    
    if (!unknown_chroms.empty) {
        throw new Exception("Unknown chromosomes found in intervals: ${unknown_chroms.unique().join(', ')}")
    }
    
    def shards = []
    chromosomes.each { chrom ->
        if (chromosomes_data.containsKey(chrom)) {
            def shard_file = file("${chrom}.bed")
            shard_file.write(chromosomes_data[chrom].join('\n') + '\n')
            shards.add([chrom, shard_file])
        }
    }
    return shards
}

process gatk_mutect_shard {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory {Math.min(56, 36+(10 * (task.attempt-1))).GB}
    cpus 12

    input:
        tuple val(sample_pair_key), val(chromosome), path(normal_bam), path(normal_bam_index), path(tumor_bam), path(tumor_bam_index)
        path(panel_of_normals)
        path(panel_of_normals_index)
        path(germline_resource)
        path(germline_resource_index)
        path(refgenome)
        path(refgenome_index)
        path(refgenome_dict)
        path(intervals)

    output:
        tuple val(sample_pair_key), val(chromosome), path("${sample_pair_key}_unfiltered_${chromosome}.vcf"), path("${sample_pair_key}_unfiltered_${chromosome}.vcf.stats"), emit: vcf_shard
        tuple val(sample_pair_key), val(chromosome), path("${sample_pair_key}_f1r2_data_${chromosome}.tar.gz"), emit: f1r2_shard

    script:
    def out_prefix = "${sample_pair_key}_${chromosome}"
    """
    # Create BED file for this chromosome
    grep -E "^${chromosome}\\s" ${intervals} > "${chromosome}.bed"
    
    mkdir -p tmp
    
    gatk Mutect2 \\
        --java-options "-Djava.io.tmpdir=tmp -Xms45G -Xmx45G" \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir tmp \\
        --input ${tumor_bam} \\
        --input ${normal_bam} \\
        -normal ${normal_bam.getSimpleName()} \\
        -L "${chromosome}.bed" \\
        -pon ${panel_of_normals} \\
        -germline-resource ${germline_resource} \\
        --reference ${refgenome} \\
        --f1r2-tar-gz "${sample_pair_key}_f1r2_data_${chromosome}.tar.gz" \\
        --output "${sample_pair_key}_unfiltered_${chromosome}.vcf"
    """
}


process gatk_merge_mutect_shards {
    conda "bioconda::gatk4=4.4.0.0"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    memory '32 GB'

    input:
        tuple val(sample_pair_key), path(vcf_shards), path(stats_shards), path(f1r2_shards)

    output:
        tuple path("${sample_pair_key}_unfiltered.vcf"), 
              path("${sample_pair_key}_unfiltered.vcf.stats"), emit: merged_vcf
        path("${sample_pair_key}_f1r2_data.shards.tar"), emit: merged_f1r2

    script:
    def vcf_names = vcf_shards.collect { it.getName() }.join('\n')
    def stats_names = stats_shards.collect { it.getName() }.join('\n')
    def f1r2_names = f1r2_shards.collect { it.getName() }.join(' ')
    """
    echo -e "${vcf_names}" > vcf_inputs.list
    echo -e "${stats_names}" > stats_inputs.list
    
    gatk MergeVcfs -I vcf_inputs.list -O "${sample_pair_key}_unfiltered.vcf"
    gatk MergeMutectStats --stats stats_inputs.list -O "${sample_pair_key}_unfiltered.vcf.stats"
    
    tar -cf "${sample_pair_key}_f1r2_data.shards.tar" ${f1r2_names}
    """
}

workflow variant_call_sharded {
  take:
    bam_pairs_w_idx
    intervals
    panel_of_normals
    panel_of_normals_index
    germline_resource
    germline_resource_index
    refgenome
    refgenome_index
    refgenome_dict

  main:
    // Create chromosome list channel
    def chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                      'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                      'chr20', 'chr21', 'chr22', 'chrX', 'chrY'] //, 'chrM']
    
    chromosome_list = Channel.fromList(chromosomes)

    // Combine: each BAM pair with each chromosome
    bam_chromosome_combinations = bam_pairs_w_idx
        .combine(chromosome_list)
        .map{ it ->
            def normal_bam = it[0]
            def normal_bai = it[1] 
            def tumor_bam = it[2]
            def tumor_bai = it[3]
            def chromosome = it[4]
            def sample_pair_key = vcf_name_from_pair(normal_bam, tumor_bam)
            
            [sample_pair_key, chromosome, normal_bam, normal_bai, tumor_bam, tumor_bai]
        }

    // Run Mutect2 shards
    mutect_shards = gatk_mutect_shard(
        bam_chromosome_combinations,
        panel_of_normals,
        panel_of_normals_index,
        germline_resource,
        germline_resource_index,
        refgenome,
        refgenome_index,
        refgenome_dict,
        intervals
    )

    // Group shards by sample pair key using groupTuple with size
    def expected_chromosome_count = 24 // chr1-chr22, chrX, chrY, chrM
    
    grouped_vcf_shards = mutect_shards.vcf_shard
        .groupTuple(by: 0, size: expected_chromosome_count)
        .map{ it ->
            def key = it[0]
            def chromosomes_list = it[1]  
            def vcfs = it[2]        
            def stats = it[3]      
            [key, vcfs, stats]
        }

    grouped_f1r2_shards = mutect_shards.f1r2_shard
        .groupTuple(by: 0, size: expected_chromosome_count)
        .map{ it ->
            def key = it[0]
            def chromosomes_list = it[1]  // list of chromosome names
            def f1r2_files = it[2]        // list of file tuples
            [key, f1r2_files]
        }

    // Combine VCF and F1R2 shards
    combined_shards = grouped_vcf_shards
        .join(grouped_f1r2_shards)
        .map{ it ->
            def key = it[0]
            def vcfs = it[1]
            def stats = it[2]
            def f1r2_files = it[3]
            [key, vcfs, stats, f1r2_files]
        }


    // Merge shards
    merged_results = gatk_merge_mutect_shards(combined_shards)

  emit:
    vcf = merged_results.merged_vcf
    f1r2_data = merged_results.merged_f1r2
}




workflow variant_call {
  take:
    bam_pairings
    sample_bams
    intervals
    panel_of_normals
    germline_resource
    refgenome
    args

  main:

    // best practices from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
    // How to Call somatic mutations using GATK4 Mutect2

    // key is XXX_T_[FFPE|BLOOD|FF]_1
    sample_bams_w_key = sample_bams.map{ it -> ["${it.getSimpleName()}", it]}
    sample_bams_idx_w_key = index_bam(sample_bams).map{ it -> ["${it.getSimpleName()}", it] }
    // join bams and bai files
    sample_bams_w_indices_w_key = sample_bams_w_key.join(sample_bams_idx_w_key)
    sample_bams_w_indices = sample_bams_w_indices_w_key.map{it -> it[1..-1]}


    bams_split = sample_bams_w_indices_w_key.branch { key, bam, bai ->
	// Match pattern: <alphanum>-<2-digit>_<N or T>_
	def normal = key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_N_.+/
	def tumor  = key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_T_.+/

	// Output channels
	normal: normal
	tumor:  tumor
    }

    // Join with normal bams: match normal_key == key
    with_normal = bam_pairings
	.combine(bams_split.normal, by: 0)  // emits [normal_key, tumor_key, normal_bam, normal_bai]
	.map { normal_key, tumor_key, normal_bam, normal_bai ->
	    [tumor_key, normal_bam, normal_bai]  // prepare for tumor join
	}

    // Join with tumor bams: match tumor_key == key
    bam_pairs_w_idx = with_normal
	.combine(bams_split.tumor, by:0)  // finds matching tumor
	.map { tumor_key, normal_bam, normal_bai, tumor_bam, tumor_bai ->
	    [normal_bam, normal_bai, tumor_bam, tumor_bai]
	}

    // join then branch the [key, [bam, bai]]

    // join then its [simplename_normal, [simplename_normal, simplename_tumor, normal bam, normal bai]]
    // map away the simplename_normal, then join on the simplename_tumor
    /*
    bam_pairs_w_idx = bam_pairings.join( sample_bams_w_indices_w_key, remainder:false).map{it -> it[1..-1]}.join(sample_bams_w_indices_w_key)
        // also map away simplename_tumor, remainder is [normal.bam, normal.bai, tumor.bam, tumor.bai]
        .map { it -> it[1..-1]}
    */

    refgenome_index = index_fasta(args.refgenome).fasta_index
    refgenome_dict = gatk_createsequencedictionary(args.refgenome).refgenome_dict

    indices = Channel.of(germline_resource, panel_of_normals)
    indices_mid = gatk_indexfeaturefile(indices).known_sites_index

    // panel_of_normals_index = gatk_indexfeaturefile(panel_of_normals).known_sites_index
    panel_of_normals_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(panel_of_normals).getSimpleName()}" }
    germline_resource_index = indices_mid.first{ it -> "${file(it).getSimpleName()}" == "${file(germline_resource).getSimpleName()}" }


def shard_mutect = true
	if(shard_mutect){
	    mut = variant_call_sharded(bam_pairs_w_idx, intervals, 
		    panel_of_normals, panel_of_normals_index,
		    germline_resource,  germline_resource_index,  
		    refgenome, refgenome_index, refgenome_dict)
	} else {
	    mut = gatk_mutect(bam_pairs_w_idx, intervals, 
		    panel_of_normals, panel_of_normals_index,
		    germline_resource,  germline_resource_index,  
		    refgenome, refgenome_index, refgenome_dict)
	}

    vcf = mut.vcf.map{vcf, stats -> vcf}
    om = gatk_learn_readorientationmodel(mut.f1r2_data).orientation_model

    // for tumors and for normals
    // todo, check that ew can use the same germline resource for pileup here

    all_pileups = gatk_getpileupsummaries(sample_bams_w_indices, intervals, germline_resource, germline_resource_index).table

    pileup_pairings = bam_pairings //.map{ normal_bam_name, tumor_bam_name -> ["${normal_bam_name}_pileup.table", "${tumor_bam_name}_pileup.table"] }

    pileups_with_key = all_pileups.map { pileup ->
	def key = removeSuffix(pileup.getName(), "_pileup.table") // e.g., TEST-24_N_FFPE_1_pileup.table -> TEST-24_N_FFPE_1
	[key, pileup]
    }

    pileup_split = pileups_with_key.branch {
	key, pileup ->
	normal: key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_N_.+/
	tumor:  key ==~ /^[a-zA-Z0-9]+-[0-9]{2}_T_.+/
    }

    // Join pairing with normal pileup
    with_normal_pileup = bam_pairings
	.combine(pileup_split.normal, by:0)  // join on normal_key == key
	.map { normal_key, tumor_key, normal_pileup ->
	    [tumor_key, normal_pileup]
	}

    // Join with tumor pileup
    paired_pileups = with_normal_pileup
	.combine(pileup_split.tumor, by:0)  // join on tumor_key == key
	.map { tumor_key, normal_pileup, tumor_pileup ->
	    [normal_pileup, tumor_pileup]
	}

    contamination_table_and_segments = gatk_calculate_contamination(paired_pileups).contamination_table_and_segments
    c_w_key = contamination_table_and_segments.map{normal_recal, tumor_recal -> 
	[
	removeSuffix(normal_recal.getSimpleName(),'_contamination'),
	normal_recal, tumor_recal
        ]}
    vcf_w_stats_w_key = mut.vcf.map{vcf, stats -> ["${removeSuffix(vcf.getSimpleName(), '_unfiltered')}", [vcf, stats]]}
    om_w_key = om.map{it -> ["${removeSuffix(it.getSimpleName(),"_f1r2_data_model")}", it]}


    ccc = vcf_w_stats_w_key.combine(c_w_key, by: 0).combine(om_w_key, by: 0)
    ccc.view()

    /*
    vcf_w_filter_data = ccc.map{ key, contamination_table, orientation_model -> 
        [contamination_table[0], contamination_table[1], orientation_model[0], orientation_model[1], orientation_model[3]
        ]} 
        */

    vcf_w_filter_data = ccc.map{ it -> [it[1][0], it[1][1], it[2], it[3], it[4]]}

    filtered_vcf = gatk_filter_calls(vcf_w_filter_data, args.refgenome, refgenome_index, refgenome_dict).vcf

  emit:
    bam_pairs_w_idx = bam_pairs_w_idx
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
