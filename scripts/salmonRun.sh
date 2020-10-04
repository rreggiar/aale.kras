#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName  directory containing trimmed *.fq.gz reads"
    echo "example $scriptName {ctrl,kras}.{1,2,3..} /path/to/{name.of.salmon.index} pipeline dateTime parameter"
    exit 1
fi

dateStamp="$3"
set -x
echo "script: $scriptName"
echo "time: $dateStamp"
set +x 

inputDir="$1"
salmonIndex="$2"
outputDir=$(basename "$salmonIndex")_${dateStamp}_out

set -x
echo "input: $inputDir"
echo "output: $outputDir"
set +x 

## activate correct env
function condaCheck() {
	# source the conda script so this shell has access
	source /public/groups/kimlab/.install_bin/anaconda3/etc/profile.d/conda.sh

	reqEnv="aale.analysis.env"
	env=$(basename "$CONDA_PREFIX")

	if [[ env != reqEnv ]]; then
		echo "switching from "$env" to "$reqEnv""
		conda activate $reqEnv
	else
		echo ""$reqEnv" is active"
	fi
}

function runSalmon() {
	inputDir="$1"
	salmonIndex="$2"
	outputDir="$3"
	outputPath="$inputDir"/"$outputDir"

	if [ ! -f "$outputPath"/quant.sf ]; then

		mkdir "$outputPath"

		trim_fwd="$inputDir"/output_forward_paired.fq.gz
		trim_rev="$inputDir"/output_reverse_paired.fq.gz

		set -x

		salmon quant \
		-i "$salmonIndex" \
		--libType A \
		-1 "$trim_fwd" \
		-2 "$trim_rev" \
		-p 8 \
		--validateMappings \
		--gcBias \
		--seqBias \
	    --recoverOrphans \
	    --rangeFactorizationBins 4 \
		--output "$outputPath" 

		set +x

	fi
}

condaCheck

runSalmon "${inputDir}" "${salmonIndex}" "${outputDir}"

