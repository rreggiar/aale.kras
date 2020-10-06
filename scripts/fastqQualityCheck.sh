#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 1 ]; then
    echo "error: usage "$scriptName" sampleDir"
    echo "example "$scriptName" /scratch/kimlab/projects/exoRNA-biomarkers-panc/data/a549_0.2MOI_24hr"
    exit 1
fi

sampleDir="$1"

dateStamp="$(bash dateStamp.sh)"

echo "cmd: "$@""

for inputDir in "$sampleDir"/*; do

	set -x 

	fastqc -t 8 output*

	exitStatus=$?
	if [ $? -ne 0 ]; then

	    echo ERROR "$scriptName" "$inputDir" returned exit status "$exitStatus"
	    continue

	fi

	set +x

done

multiqc "$sampleDir" --ingore *.fastq.gz --filename multiqc."${dateStamp}" -outdir "$sampleDir"