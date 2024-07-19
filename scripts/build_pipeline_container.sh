#!/bin/bash

containerfile=pipeline_container.containerfile
tag="${containerfile%%.containerfile}"
podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" .
