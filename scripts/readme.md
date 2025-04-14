this folder contains scripts to build a container that includes all dependencies for the pipeline.

it uses the build_deps.sh and install_deps.sh scripts.

those could otherwise be used to manually install all dependencies locally.

the tools are built from their github sources, to ease development and bevcause they're not all packaged yet.

the pipeline_task container can be built with:

```
TMPDIR=./tmp ./gen_containerfile.sh
```

the reference data can be downloaded with the download_igenomes.py script

the script only depends on boto3


the reference files can also be obtained with: https://ewels.github.io/AWS-iGenomes/

the hashes can be verified in parallel with:

```
pip install boto3
python download_igenomes.py PATH_WHERE_TO_STORE_REFERENCES
cp igenomes_sha256sum.txt PATH_WHERE_TO_STORE_REFERENCES
cd PATH_WHERE_TO_STORE_REFERENCES
cat igenomes_sha256sum.txt | parallel --pipe -N1 sha256sum -c
```
