#!/bin/bash
echo "the testdata will be downloaded to TESTDATA_DIR"

if [ -d "$TESTDATA_DIR" ]; then
	#wget https://www.ncbi.nlm.nih.gov/sra/SRP162370 --output-document=SRR7890827.sra

	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890824/SRR7890824 --output-document=$TESTDATA_DIR/SRR7890824.sra
	wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890827/SRR7890827 --output-document=$TESTDATA_DIR/SRR7890827.sra
else
	echo "error, the env-var TESTDATA_DIR must be set to an existing directory"
fi
