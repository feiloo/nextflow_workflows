nextflow.enable.dsl=2

process sambamba_markdup {
    conda "sambamba=1.0.1"
    container 'quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_0'

    cpus { Math.max(1, Math.round(26 * (1 - ((1/4)*(task.attempt-1))))) }
    errorStrategy 'retry'
    maxRetries 4

    memory {(200 * (1.3 ** (task.attempt-1))).GB}

    input:
        path(bam)

    output:
      path("out/${bam}"), emit: marked_bams

    script:
    def compression_level = 1

    """
    mkdir out
    mkdir tmp
    sambamba markdup --tmpdir tmp -l ${compression_level} -t ${task.cpus} ${bam} out/${bam}
    """
}
