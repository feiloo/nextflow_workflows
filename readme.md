# Readme
env-vars:

export NXF_HOME='$HOME/.nextflow'

export NXF_ASSETS='$NXF_HOME/assets'

export NXF_OPTS='-Dlog=/path/nextflow_logs'

export NXF_TEMP=/path/nextflow_temp/

export NXF_LOG_FILE='/path/nextflow_logs'

export NXF_PLUGINS_DIR='/path/nextflow_plugins'

export NEXTFLOW_MODULES="$(pwd)/modules"

### build with meson

exporting the right NEXTFLOW_MODULES env var is sufficient for running these nextflow modules.

but the meson buildsystem enables additional steps for building, testing and packaging.

#### create a builddir

meson uses the advantageous out of source build directories
create one with:

meson setup --wipe /path_to_builddir/

#### compile with meson

cd /path_to_builddir/
meson compile

#### run tests with meson

cd /path_to_builddir
meson test
