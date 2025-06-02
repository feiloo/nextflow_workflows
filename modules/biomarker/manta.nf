nextflow.enable.dsl=2

process manta_somatic {
    
    tag "${sample_id}"
    debug true
            
    conda "${moduleDir}/environment.yml"
  
    publishDir(path: "${outdir}/${sample_id}/manta_results", mode: "copy")    
     
    input:
        tuple val(sample_id), path(normal_bam), path(normal_bai), path(tumor_bam) , path(tumor_bai)
        tuple val(id1), val(fasta)
        tuple val(id2), val(fai)
        val(call_region_grch38_bed)
        val(manta_config_file)
        
    output:
        tuple val(sample_id), path("${sample_id}.candidate_small_indels.vcf.gz")     , emit: candidate_small_indels_vcf
        tuple val(sample_id), path("${sample_id}.candidate_small_indels.vcf.gz.tbi") , emit: candidate_small_indels_vcf_tbi
        tuple val(sample_id), path("${sample_id}.candidate_sv.vcf.gz")               , emit: candidate_sv_vcf
        tuple val(sample_id), path("${sample_id}.candidate_sv.vcf.gz.tbi")           , emit: candidate_sv_vcf_tbi
        tuple val(sample_id), path("${sample_id}.diploid_sv.vcf.gz")                 , emit: diploid_sv_vcf
        tuple val(sample_id), path("${sample_id}.diploid_sv.vcf.gz.tbi")             , emit: diploid_sv_vcf_tbi
        tuple val(sample_id), path("${sample_id}.somatic_sv.vcf.gz")                 , emit: somatic_sv_vcf
        tuple val(sample_id), path("${sample_id}.somatic_sv.vcf.gz.tbi")             , emit: somatic_sv_vcf_tbi
        tuple val(sample_id), path("${sample_id}.alignmentStatsSummary.txt")         , emit: alignmentStatsSummary 
        tuple val(sample_id), path("${sample_id}.svCandidateGenerationStats.tsv")    , emit: svCandidateGenerationStats_tsv
        tuple val(sample_id), path("${sample_id}.svCandidateGenerationStats.xml")    , emit: svCandidateGenerationStats_xml 
        tuple val(sample_id), path("${sample_id}.svLocusGraphStats.tsv")             , emit: svLocusGraphStats
        tuple val(sample_id), path("${sample_id}.manta_workflow.log")                , emit: manta_workflow_log
                
    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    """
    configManta.py --normalBam ${normal_bam} \\
                   --tumorBam ${tumor_bam} \\
                   --referenceFasta ${fasta} \\
                   --callRegions ${call_region_grch38_bed} \\
                   --config ${manta_config_file} \\
                   --runDir manta

    python manta/runWorkflow.py -m $args -j $args2 > ${sample_id}.manta_workflow.log 2>&1

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${sample_id}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${sample_id}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${sample_id}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${sample_id}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        ${sample_id}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        ${sample_id}.diploid_sv.vcf.gz.tbi
    mv manta/results/variants/somaticSV.vcf.gz \\
        ${sample_id}.somatic_sv.vcf.gz
    mv manta/results/variants/somaticSV.vcf.gz.tbi \\
        ${sample_id}.somatic_sv.vcf.gz.tbi
    mv manta/results/stats/alignmentStatsSummary.txt \\
        ${sample_id}.alignmentStatsSummary.txt
    mv manta/results/stats/svCandidateGenerationStats.tsv \\
        ${sample_id}.svCandidateGenerationStats.tsv
    mv manta/results/stats/svCandidateGenerationStats.xml \\
        ${sample_id}.svCandidateGenerationStats.xml
    mv manta/results/stats/svLocusGraphStats.tsv \\
        ${sample_id}.svLocusGraphStats.tsv       
    """
}

csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, file(row.normal_bam), file(row.tumor_bam)) }

//csv_ch.view()
ref_fasta_ch = Channel.value(params.ref_fasta)
ref_fasta_fai_ch = Channel.value(params.ref_fasta_fai)
call_region_grch38_bed_ch = Channel.value(params.call_region_grch38_bed)
manta_config_file_ch = Channel.value(params.manta_config_file)
outdir = params.outdir
moduleDir = params.moduleDir

workflow manta {

    take:
       csv_ch
       ref_fasta_ch
       ref_fasta_fai_ch
       call_region_grch38_bed_ch
       manta_config_file_ch     

    main:
       manta_bam_index(csv_ch)
       
       fasta_ch = ref_fasta_ch.map { file -> 
                                def id = "id1"
                                def file_path = file

                                tuple(id, file_path)
                         }
       fai_ch = ref_fasta_fai_ch.map { file_fai ->
                                def id_fai = "id2"
                                def file_fai_path = file_fai

                                tuple(id_fai, file_fai_path)
                         }
       //fasta_fai_ch = fasta_ch.combine(fai_ch, by: 0)
       //manta_somatic(manta_bam_index.out,params.ref_fasta,params.ref_fasta_fai,call_region_grch38_bed_ch)
       manta_somatic(manta_bam_index.out,fasta_ch,fai_ch,call_region_grch38_bed_ch,manta_config_file_ch)      
 }

workflow {
    manta(csv_ch,ref_fasta_ch,ref_fasta_fai_ch,call_region_grch38_bed_ch,manta_config_file_ch)
    //manta(csv_ch,call_region_grch38_bed_ch)
}
