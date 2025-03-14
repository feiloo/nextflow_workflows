process fastp {
    conda "bioconda::fastp=0.24.0"
    container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    cache 'lenient'

    //cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    cpus 8
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(sample_id), path(read1), path(read2), val(output_file_prefix)

    output:
    tuple val(sample_id), path("${output_file_prefix}_fastp.html"), emit: html
    tuple val(sample_id), path("${output_file_prefix}_fastp.json"), emit: json
    tuple val(sample_id), path("out/${read1.getSimpleName()}.fq.gz"), path("out/${read2.getSimpleName()}.fq.gz"), emit: preprocessed_reads

    script:

    def gz_compressionlevel = 9 // out of 1 to 9

    """
    mkdir -p out

    fastp \\
	--in1 ${read1} \\
	--in2 ${read2} \\
	--out1  out/${read1.getSimpleName()}.fq.gz \\
	--out2  out/${read2.getSimpleName()}.fq.gz \\
	-z ${gz_compressionlevel} \\
	--thread ${task.cpus} \\
	--json ${output_file_prefix}_fastp.json \\
	--html ${output_file_prefix}_fastp.html \\
	2> ${output_file_prefix}.fastp.log
    """
}
