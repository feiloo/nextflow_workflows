process fastp {
    conda "bioconda::fastp=1.0.1"
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'

    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    cache 'lenient'

    //cpus { Math.max(1, Math.round(Runtime.runtime.availableProcessors() * (1 - ((1/4)*(task.attempt-1))))) }
    cpus 12
    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4

    input:
    tuple val(sample_id), path(read1), path(read2), val(output_file_prefix), val(expected_sha256sum1), val(expected_sha256sum2)

    output:
    tuple val(sample_id), path("${output_file_prefix}_fastp.html"), emit: html
    tuple val(sample_id), path("${output_file_prefix}_fastp.json"), emit: json
    tuple val(sample_id), path("out/${read1.getSimpleName()}.fq.gz"), path("out/${read2.getSimpleName()}.fq.gz"), emit: preprocessed_reads
    path("${output_file_prefix}_check_process_results.txt"), emit: integrity_check

    script:

    def gz_compressionlevel = 9 // out of 1 to 9

    """
    set -euo pipefail

    mkdir -p out

    fastp_preprocessing() {
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
	}


    check_integrity() {
	    NAME=\$1
	    NAMEOUT="\$NAME"_out
	    expected_sha256sum=\$2

	    mkdir "\$NAME"_out
	    sha256sum \$NAME | awk '{print \$1}' > \$NAMEOUT/sha256sum_result &
	    SHA256SUM_PID=\$!

	    #md5sum \$NAME > \$NAMEOUT/md5sum_result

	    gzip -t --keep \$NAME > \$NAMEOUT/gzip_result &
	    GZIP_PID=\$!
	    bgzip -t -@ 8 --keep \$NAME > \$NAMEOUT/bgzip_result
	    BGZIP_PID=\$!


	    wait \$SHA256SUM_PID

	    SHA256SUM=\$(cat \$NAMEOUT/sha256sum_result)
	    if [ "\$SHA256SUM" != "\${expected_sha256sum}" ]; then
		echo "error, sha256sum of \$NAME doesnt match: \$SHA256SUM  \$expected_sha256sum"
		exit 250

	    fi 
	    #MD5SUM_PID=\$!
	    wait \$GZIP_PID

	    if grep -qw 'trailing' gzip_result && grep -qw 'garbage' gzip_result; then
		echo "error, found trailing garbage in \$NAME according to gzip. the file is possibly truncated"
		exit 251
	    fi

	    # todo add verficiation of bgzip too here
	    wait \$BGZIP_PID

	    echo -e "check process results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/sha256sum_result >> \$NAMEOUT/check_process_results
	    echo -e "\ngzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/gzip_result >> \$NAMEOUT/check_process_results
	    echo -e "\nbgzip check results:\n" >> \$NAMEOUT/check_process_results.txt
	    cat \$NAMEOUT/bgzip_result >> \$NAMEOUT/check_process_results

    }

    fastp_preprocessing &
    FASTP_PID=\$!

    check_integrity "${read1}" "${expected_sha256sum1}" &
    CHECK1=\$!
    check_integrity "${read2}" "${expected_sha256sum2}" &
    CHECK2=\$!

    wait \$CHECK1
    wait \$CHECK2
    wait \$FASTP_PID

    cat ${read1}_out/check_process_results ${read2}_out/check_process_results > ${output_file_prefix}_check_process_results.txt

    """
}

/*
// prototype for separate integrity check process

process integrity_check {
    conda "bioconda::bcftools=1.17"
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"

    // check gzip and bgzip archive intergrity

    cache 'lenient'
    cpu 10
    memory '1 GB'

    // set to "ignore" to filter out broken files
    // or to "terminate" to kill the whole run
    errorStrategy "ignore"

    input: 
    tuple path(fastq), val(expected_sha256sum)

    output:
    tuple path("${fastq}_check_process_results")

    script:
    """
    # could optimize to read data only once
    # cat ${fastq} | tee >(gzip -t --keep > gzip_result) | tee >(bgzip -@ 2 -t > bgzip_result) | tee >(sha256sum > sample_sha256sum) | md5sum > sample_md5sum

    # but filesystem cache could be good enough
    sha256sum ${fastq} | awk '{print \$1}' > sha256sum_result
    SHA256SUM_PID=\$!

    #md5sum ${fastq} > md5sum

    gzip -t --keep ${fastq} > gzip_result &
    GZIP_PID=\$!
    bgzip -t -@ 8 --keep ${fastq} > bgzip_result
    BGZIP_PID=\$!

    wait \$SHA256SUM_PID

    SHA256SUM=\$(cat sha256sum_result)
    if [ "\$SHA256SUM" != "${expected_sha256sum}" ]; then
        echo "error, sha256sum doesnt match"
        exit 1

    if 
    #MD5SUM_PID=\$!
    wait \$GZIP_PID

    if grep -qw 'trailing' gzip_result && grep -qw 'garbage' gzip_result; then
        echo "error, found trailing garbage according to gzip. the file is possibly truncated"
	exit 2
    fi

    # todo add verficiation of bgzip too here
    wait \$BGZIP_PID

    echo -e "check process results:\n" >> check_process_results.txt
    cat sha256sum_result >> check_process_results
    echo -e "\ngzip check results:\n" >> check_process_results.txt
    cat gzip_result >> check_process_results
    echo -e "\nbgzip check results:\n" >> check_process_results.txt
    cat bgzip_result >> check_process_results

    mv check_process_results ${fastq}_check_process_results

    #wait $MD5SUM_PID
    """

}
*/
