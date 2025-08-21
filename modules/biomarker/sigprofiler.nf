nextflow.enable.dsl=2
   
process extractor {

    tag "${sample_id}"
    debug true

    publishDir(path: "${outdir}/sigprofiler_results/", mode: "copy")

    input:
    tuple val(sample_id), path(vcf)

    output:
    path("*")

    script:
    """
    SigProfilerExtractor sigprofilerextractor vcf ${sample_id} ${vcf} --reference_genome GRCh38  --opportunity_genome GRCh38 > ${sample_id}.run.log
    """
}

csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sample_id, file(row.vcf_path)) }
outdir = params.outdir

workflow sigprofiler {

    take:
       csv_ch
              
   main:
       extractor(csv_ch)
}

workflow {
    sigprofiler(csv_ch)
}
