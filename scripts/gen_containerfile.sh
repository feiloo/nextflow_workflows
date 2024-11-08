#!/bin/bash
set -euo pipefail

# generate the sub-scripts that make up the steps of the containerfile and build the container
# this allows for layered caching of the container build, which helps because the build is slow

temp_dir=$(mktemp -d ${TMPDIR:-/tmp}/nextflow_generated_scripts.XXXXXX)
echo $temp_dir

cp build_deps.sh $temp_dir
cp install_deps.sh $temp_dir

# copy the template containerfile into output
cp pipeline_task.containerfile $temp_dir

# split the bash script into buildsteps-scripts
pushd $temp_dir
cat build_deps.sh | csplit - /###buildstep/ {*}
chmod ug+rx *.sh
chmod ug+rx xx*

echo "RUN mkdir gen_containerfiles" >> pipeline_task.containerfile

# generate containerfile instructions that copy and run the scripts within the container
ls xx* | xargs -i echo -e "COPY {} /root/gen_containerfiles \nRUN /root/gen_containerfiles/{}" >> pipeline_task.containerfile

echo -e "COPY install_deps.sh /root/\n" >> pipeline_task.containerfile
echo -e "RUN chmod u+x install_deps.sh && ./install_deps.sh /root/\n" >> pipeline_task.containerfile

pushd $temp_dir

containerfile=pipeline_task.containerfile
tag="${containerfile%%.containerfile}"
echo "$containerfile" container building from $temp_dir
podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" .

popd
