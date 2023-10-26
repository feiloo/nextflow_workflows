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

```
cd /path_to_builddir
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
