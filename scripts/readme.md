this folder contains scripts to build a container that includes all dependencies for the pipeline.

it uses the build_deps.sh and install_deps.sh scripts.

those could otherwise be used to manually install all dependencies locally.

the tools are built from their github sources, to ease development and bevcause they're not all packaged yet.

the pipeline_task container can be built with:

```
TMPDIR=./tmp ./gen_containerfile.sh
```
