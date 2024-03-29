process fastp {
    conda "bioconda::fastp=0.23.4"
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    memory "6 GB"

    input:
    tuple val(sample_id), path(read1), path(read2), val(output_file_prefix)

    output:
    tuple val(sample_id), path("${output_file_prefix}_fastp.html"), emit: html
    tuple val(sample_id), path("${output_file_prefix}_fastp.json"), emit: json
    tuple val(sample_id), path("out/${read1.getSimpleName()}.fq.gz"), path("out/${read2.getSimpleName()}.fq.gz"), emit: preprocessed_reads

    script:

    def gz_compressionlevel = 1 // out of 1 to 9

    n_cpus = Runtime.runtime.availableProcessors()
    """
    mkdir out

    fastp \\
	--in1 ${read1} \\
	--in2 ${read2} \\
	--out1  out/${read1.getSimpleName()}.fq.gz \\
	--out2  out/${read2.getSimpleName()}.fq.gz \\
	-z ${gz_compressionlevel} \\
	--thread $n_cpus \\
	--json ${output_file_prefix}_fastp.json \\
	--html ${output_file_prefix}_fastp.html \\
	2> ${output_file_prefix}.fastp.log
		
    """
}
