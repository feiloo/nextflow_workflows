process fastqc {
    conda "bioconda::fastqc=0.11.9"
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    memory "60 GB"

    cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${read.getSimpleName()}_fastqc.html"), emit: html
    tuple val(sample_id), path("${read.getSimpleName()}_fastqc.zip"), emit: zip

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    fastqc \\
        --threads ${task.cpus} \\
	${read}
        
    """
}
