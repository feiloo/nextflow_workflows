#!/bin/bash
podman run -it --rm \
	-v /var/run/docker.sock:/var/run/docker.sock \
	-v ${TMPDIR:-/tmp}:${TMPDIR:-/tmp} \
	nextflow
