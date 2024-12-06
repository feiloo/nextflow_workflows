# Nextflow modules for NGS

this repo contains multiple nextflow workflows for somatic NGS analysis.

we focus on modularity and simplicity.

the implementation is inspired by various nf-core modules and the nf-core/sarek pipeline.

## Setup

there are multiple ways of running the workflow. 
the simplest way is to run (portable) from a single folder, but WARNING this pipeline requires a few terabytes of free storage.

for now we recommend running on a single host with conda.

Dependencies:

* nextflow
* conda

```
git clone --recursive https://github.com/feiloo/nextflow_workflows.git
```

#### Set required environment variables

```
mkdir reference nextflow_calldir nextflow_workdir nextflow_outputdir cache
export NEXTFLOW_MODULES="$(pwd)/nextflow_workflows/modules"
export NGS_REFERENCE_DIR="$(pwd)/reference"
export NEXTFLOW_CALLDIR="$(pwd)/nextflow_calldir"
export NEXTFLOW_WORKDIR_CUSTOM="$(pwd)/nextflow_workdir'
export NEXTFLOW_OUTPUTDIR_CUSTOM="$(pwd)/nextflow_outputdir"
export NEXTFLOW_STOREDIR="$(pwd)/cache"
```

#### Download required reference data

```
# requires boto3
python scripts/download_data.py reference $NGS_REFERENCE_DIR
```

## usage

the pipeline has multiple variations, the main variation is "align_interpret"

it requires a samplesheet.csv like this:

```
sample_id,normal_read1,normal_read2,tumor_read1,tumor_read2
SAMPLENAME-24,/data/testdata/SAMPLENAME-24_N_1.fq.gz,/data/testdata/SAMPLENAME-24_N_2.fq.gz,/data/testdata/SAMPLENAME-24_T_1.fq.gz,/data/testdata/SAMPLENAME-24_T_2.fq.gz
```

now start the pipeline.

```
nextflow run $NEXTFLOW_MODULES/ukb_main_workflow/main.nf \
    -profile standard \
    --workflow_variation align_interpret \
    --samplesheet test_sequence_alignment_samplesheet.csv
```

## configuration

see the environment variables and nextflow-configs like modules/ukb_main_workflow/user.config


## development

for development, you can set the NEXTFLOW_MODULES variable to this repos modules directory:
```
export NEXTFLOW_MODULES="$(pwd)/modules"
```

development dependencies:

* meson

### building

the meson buildsystem enables additional steps for building, testing and packaging.

#### create a builddir

meson uses out of source build directories, which has various advantages.

create one with:

```
meson setup --wipe -D test_datadir=path_to_datadir/ -D external_deps=path_to_external_deps/ path_to_builddir/
```

#### compile with meson

```
cd /path_to_builddir/
meson compile
```

#### run tests with meson

test setup: todo

now configure the tests with them:

```
cd /path_to_builddir
meson configure -Dtest_datadir=/path_to_test_samplesheet_dir
```

run tests

```
meson test
```

#### install tools natively

instead of conda or containers, the required tools within the pipeline can be installed locally, see:

```
scripts/build_deps.sh
scripts/install_deps.sh
```

#### use merged container

another way is to use a single container for all processes or for the whole pipeline, see:

```
scripts/gen_container.sh
```


## Other notes and hints

### useful nextflow environment variables:

```
export NXF_HOME='$HOME/.nextflow'
export NXF_ASSETS='$NXF_HOME/assets'
export NXF_OPTS='-Dlog=/path/nextflow_logs'
export NXF_TEMP=/path/nextflow_temp/
export NXF_LOG_FILE='/path/nextflow_logs'
export NXF_PLUGINS_DIR='/path/nextflow_plugins'
export NEXTFLOW_MODULES="$(pwd)/modules"
```

### our paralellism strategy:

*subject to change*

use nextflow processes (dataflow) to use multinode+multiprocessing parallelise concurrently across samples/files

we work to enable multinode, but for now we just use 1 node

for now, we avoid splitting files for more parallelism

instead we use tool-level threading

so we run multiple nf-processes

each nf-process runs a single tool process 

each tool process runs the tool with the tools parallelism option


with the current filesize to count ratio: 1GB to 5 GB per file and about 5-50 files per workflow on ~100 cores this works pretty well

and most importantly it is quite simple to reason about
