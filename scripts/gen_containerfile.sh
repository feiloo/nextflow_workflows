#!/bin/bash
set -euo pipefail

# generate the sub-scripts that make up the steps of the containerfile and build the container
# this allows for layered caching of the container build, which helps because the build is slow

temp_dir=$(mktemp -d ${TMPDIR:-/tmp}/nextflow_generated_scripts.XXXXXX)
echo $temp_dir

#mkdir -p gen_containerfiles
rm -rf $temp_dir/*

cp build_deps.sh $temp_dir
cp install_deps.sh $temp_dir
mkdir -p $temp_dir/gen_containerfiles

# split the bash script into buildsteps-scripts
pushd $temp_dir/gen_containerfiles
cat ../build_deps.sh | csplit - /###buildstep/ {*}
popd

chmod u+x gen_containerfiles/x*

# copy the template containerfile into output
cp pipeline_container.containerfile $temp_dir/gen_containerfiles

echo "RUN mkdir gen_containerfiles" >> $temp_dir/gen_containerfiles/pipeline_container.containerfile

# generate containerfile instructions that copy and run the scripts within the container
ls gen_containerfiles/x* | xargs -i echo -e "COPY {} /root/gen_containerfiles \nRUN ./{}" >> $temp_dir/gen_containerfiles/pipeline_container.containerfile

pushd $temp_dir/
mkdir -p ${TMPDIR:-/tmp}
chmod ug+rx gen_containerfiles/*

echo -e "COPY install_deps.sh /root/\n" >> $temp_dir/gen_containerfiles/pipeline_container.containerfile
echo -e "RUN chmod u+x install_deps.sh && ./install_deps.sh /root/\n" >> $temp_dir/gen_containerfiles/pipeline_container.containerfile

containerfile=gen_containerfiles/pipeline_container.containerfile
tag="${containerfile%%.containerfile}"
echo building from $temp_dir
podman build --ulimit nofile=65535:65535 --tag="$tag" --file "$containerfile" .

popd
