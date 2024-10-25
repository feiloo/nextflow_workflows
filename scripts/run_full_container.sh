podman run -ti --rm \
	--mount type=bind,src=/data/reference/,destination=/usr/share/nextflow_pipeline_reference,ro=True \
	--mount type=bind,src=/data2/bwa_indices/,destination=/var/cache/nextflow_pipeline_storedir,ro=True \
	pipeline_cont
