#!/bin/bash

# aedavids@ucsc.edu
## download tools and install
# - [SRA toolkit home page](https://ncbi.github.io/sra-tools/)
# - [download tool for centos](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz)
# - [how to use sra tools](https://edwards.sdsu.edu/research/fastq-dump/)


scriptName=`basename $0`
if [ $# -lt 1 ]; then
    echo "error: usage $scriptName  list of SRA tool sample ids"
    echo "example $scriptName  SRR10080507 SRR10080508"
    exit 1
fi

set -x
sraDir=/public/groups/kimlab/NCBI/sra
listOfSRA_IDS=${@:1}

function runner() {
	listOfSampleId="$1"
	for i in $listOfSampleId
	do
	    set -x
	    # --gzip
	    # fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-e --clip ${i}

	    # download sra file
	    prefetch $i

	    # convert to fastq.gz file format
	    fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-e --clip ${sraDir}/${i}.sra
	    set +x
	done
}


runner "${listOfSRA_IDS}"

