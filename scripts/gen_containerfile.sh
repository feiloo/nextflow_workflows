#!/bin/bash
set -euo pipefail

repodir=$(pwd)/..

# generate the sub-scripts that make up the steps of the containerfile and build the container
# this allows for layered caching of the container build, which helps because the build is slow

temp_dir=$(realpath $(mktemp -d ${TMPDIR:-/tmp}/nextflow_generated_scripts.XXXXXX))
echo using tempdir: $temp_dir

# temp_dir_abs=${"$(pwd)/$temp_dir":-/tmp}

# if [ -n "$TMPDIR" ]; then
# 	temp_dir_abs=$(pwd)/$temp_dir
# else
# 	temp_dir_abs="/tmp"
# fi

cp -v build_deps.sh $temp_dir
cp -v install_deps.sh $temp_dir

# copy the template containerfile into output
cp -v pipeline_task.containerfile $temp_dir
cp -v pipeline.containerfile $temp_dir

# split the bash script into buildsteps-scripts
echo splitting and copying buildsteps scripts into tmdir
pushd $temp_dir
cat build_deps.sh | csplit - /###buildstep/ {*}
chmod ug+rx *.sh
chmod ug+rx xx*

echo generating containerfiles

echo "RUN mkdir gen_containerfiles" >> pipeline_task.containerfile

# generate containerfile instructions that copy and run the scripts within the container
ls xx* | xargs -i echo -e "COPY {} /root/gen_containerfiles \nRUN /root/gen_containerfiles/{}" >> pipeline_task.containerfile

echo -e "COPY install_deps.sh /root/\n" >> pipeline_task.containerfile
echo -e "RUN chmod u+x install_deps.sh && ./install_deps.sh /root/\n" >> pipeline_task.containerfile

echo -e "RUN chmod u+x install_deps.sh && ./install_deps.sh /root/\n" >> pipeline_task.containerfile

container_version="0.0.1"

# build the pipeline task container
containerfile=pipeline_task.containerfile
tag="${containerfile%%.containerfile}:${container_version}"
echo "$containerfile" container building from $temp_dir
TMPDIR=$temp_dir podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" .

# replace the first line in the containerfile to use the exact version of pipeline task here
sed -i "1s/.*/FROM pipeline_task:${container_version}/" pipeline.containerfile

# build the full pipeline container
containerfile=pipeline.containerfile
tag="${containerfile%%.containerfile}:${container_version}"
echo "$containerfile" container building from "$repodir"
TMPDIR=$temp_dir podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" "$repodir"

popd
