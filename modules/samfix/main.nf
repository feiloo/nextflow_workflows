process samfix {
    conda "bioconda::samtools=1.17"
    //container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    //container 'localhost/samfix'
    container 'localhost/samfix-nuitka'

    input:
    tuple path(bamfile), path(bamfile_index)

    output:
    path("out/${bamfile}"), emit: bam

    script:
    n_cpus = Runtime.runtime.availableProcessors()
    """
    mkdir -p out
    #samtools view -@ $n_cpus -bS - out/${bamfile}
    #samtools view -@ $n_cpus -h ${bamfile} | cat - > /dev/null
    #samtools view -@ $n_cpus -h ${bamfile} | python3.11 /usr/local/lib/samfix.py > out.sam
    #samtools view -@ $n_cpus -h ${bamfile} | python3.11 /usr/local/lib/samfix.py | samtools view -@ $n_cpus -o out/${bamfile} -bS -
    samtools view -@ $n_cpus -h ${bamfile} | samfix | samtools view -@ $n_cpus -o out/${bamfile} -bS -
    """
}

