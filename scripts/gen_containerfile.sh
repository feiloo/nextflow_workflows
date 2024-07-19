#!/bin/bash
set -euo pipefail

mkdir -p gen_containerfiles
rm -rf gen_containerfiles/*

# split the bash script into buildsteps-scripts
pushd gen_containerfiles
cat ../build_deps.sh | csplit - /###buildstep/ {*}
popd

chmod u+x gen_containerfiles/x*

# copy the template containerfile into output
cp pipeline_container.containerfile gen_containerfiles/

echo "RUN mkdir gen_containerfiles" >> gen_containerfiles/pipeline_container.containerfile

# generate containerfile instructions that copy and run the scripts within the container
ls gen_containerfiles/x* | xargs -i echo -e "COPY {} /root/gen_containerfiles \nRUN ./{}" >> gen_containerfiles/pipeline_container.containerfile

containerfile=gen_containerfiles/pipeline_container.containerfile
tag="${containerfile%%.containerfile}"
podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" .
