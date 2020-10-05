#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

## !!install salmon via conda!!
## $ conda config --add channels conda-forge
## $ conda install -c bioconda salmon
## $ salmon --version # double check you're at the current version

# run salmon on trimmed fastq files to quantify abundance of transcripts against provided 
# txome index. Resulting quant.sf files are used in DESeq

scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName inputDir salmonIndex dateStamp"
    echo "example $scriptName {ctrl,kras}.{1,2,3..} /path/to/{name.of.salmon.index} pipeline_dateStamp"
    exit 1
fi

dateStamp="$3"
set -x
echo "script: $scriptName"
echo "salmon version:" $(salmon --version)
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
# this will only run if you happen to be in the wrong env
# I think using directories as envNames could make this moot
# OR allow us to rely on "$(basename $PWD)" as $reqENV which would make the check
# super portable
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

