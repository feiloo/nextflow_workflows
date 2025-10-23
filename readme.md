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
export NEXTFLOW_WORKDIR_CUSTOM="$(pwd)/nextflow_workdir"
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

this version requires fastq input files with naming that matches this example:

```
X123-25_N_FFPE_1.fq.gz
X123-25_N_FFPE_2.fq.gz
X123-25_T_FFPE_1.fq.gz
X123-25_T_FFPE_2.fq.gz
```

where X123-YY is the samplename and 25 is the year.
N or T stands for normal or tumor.
FFPE or FF or BLOOD stands for formalin fixed paraffin embedded, fresh frozen, or blood source material respectively.
1 or 2 are the read direction and .fq.gz stands for bgzf/bgzip compressed fastq files.

it also requires a samplesheet.csv like this:

```
sample_id,normal_modality,tumor_modality
X123-25,FFPE,FFPE
X123-25,BLOOD,FFPE
```

a md5sum.txt that includes the hashes for the input

the pipeline takes the fastq samplesheet.csv and md5sum.txt files from the --input_dir.

now start the pipeline.

```
cd $NEXTFLOW_CALLDIR && nextflow \
        -log nextflow.log \
        run $NEXTFLOW_MODULES/ukb_main_workflow/ \
        --workflow_variation align_interpret \
        --samplesheet samplesheet.csv \
        --hash_db md5sum.txt \
        -c $NEXTFLOW_MODULES/ukb_main_workflow/user.config \
        --library_type wgs \
        --input_dir $NEXTFLOW_CALLDIR \
        -profile standard \
        -offline \
        -resume \
        --tag routine_establish,wgs
```

## configuration

see the environment variables and nextflow-configs like `modules/ukb_main_workflow/user.config`


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
