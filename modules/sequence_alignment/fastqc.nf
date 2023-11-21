process fastqc {
    conda "bioconda::fastqc=0.11.9"
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    memory "36 GB"

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${read.getSimpleName()}_fastqc.html"), emit: html
    tuple val(sample_id), path("${read.getSimpleName()}_fastqc.zip"), emit: zip

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    fastqc \\
        --threads $n_cpus \\
	${read}
        
    """
}
