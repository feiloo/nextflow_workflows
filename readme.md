# Nextflow modules for NGS

this is vey much a work in progress!

this repo contains nextflow modules which can be and are composed into NGS analysis pipelines.

we focus on modularity and simplicity.

the implementation is inspired by various nf-core modules and the nf-core/sarek pipeline.


### NEXTFLOW_MODULES environment variable

the modules in this repo include eachother from the directory in the environment variable NEXTFLOW_MODULES

modules should be installed into a sensible linux filesystem standard directory.

therefore we recommend using:
```
export NEXTFLOW_MODULES=/usr/local/lib/nextflow
```

if you want to create a portable/static workflow, you can simply copy all required modules into a single directory and set the variable.

### build with meson

for development, you can set the NEXTFLOW_MODULES variable to this repos modules directory.

the meson buildsystem enables additional steps for building, testing and packaging.

#### create a builddir

meson uses out of source build directories, which has various advantages.
create one with:

```
meson setup --wipe /path_to_builddir/
```

#### compile with meson

```
cd /path_to_builddir/
meson compile
```

#### run tests with meson

configure the tests:
```
cd /path_to_builddir
meson configure -Dtest_samplesheets=/path_to_test_samplesheet_dir -Dtest_configs=/path_to_test_configs_dir
```
within those folders create a test_main.config and a samplesheet for each workflow_variation, e.q. test_clc_samplesheet.csv

run tests
```
meson test
```

### usefull env-vars:

export NXF_HOME='$HOME/.nextflow'

export NXF_ASSETS='$NXF_HOME/assets'

export NXF_OPTS='-Dlog=/path/nextflow_logs'

export NXF_TEMP=/path/nextflow_temp/

export NXF_LOG_FILE='/path/nextflow_logs'

export NXF_PLUGINS_DIR='/path/nextflow_plugins'

export NEXTFLOW_MODULES="$(pwd)/modules"

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
