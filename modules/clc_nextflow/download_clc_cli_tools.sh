#!/bin/bash
set -euo pipefail
OUTPUT_PATH=$1
#INSTALL_SCRIPT='CLCServerCommandLineTools_24_0_1_64.sh'
INSTALL_SCRIPT=$2

if [ ! -f third_party/$INSTALL_SCRIPT ]; then
	wget https://download.clcbio.com/CLCServerCommandLineTools/24.0.1/CLCServerCommandLineTools_24_0_1_64.sh -O $OUTPUT_PATH/$INSTALL_SCRIPT
else
	echo "script file already exists, skipping download"
fi
