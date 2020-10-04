#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

scriptName=$(basename $0)
if [ $# -lt 1 ]; then
    echo "error: usage $scriptName  directory containing trimmed *.fq.gz reads"
    echo "example $scriptName  /path/to/{ctrl,kras}.{1,2,3..}/salmon.out/"
    exit 1
fi

dateStamp=$(dateStamp.sh)
echo "time: $dateStamp"

## activate correct env
# source the conda script so this shell has access
function condaCheck() {
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
	if [ ! -f "$final_output_path"/quant.sf ]; then

		mkdir "$final_output_path"

		trim_fwd="$sample"/*output_forward_paired.fq.gz
		trim_rev="$sample"/*output_reverse_paired.fq.gz

		salmon quant \
		-i "$salmon_index" \
		--libType A \
		-1 "$trim_fwd" \
		-2 "$trim_rev" \
		-p 8 \
		--validateMappings \
		--gcBias \
		--seqBias \
	    --recoverOrphans \
	    --rangeFactorizationBins 4 \
		--output "$final_output_path" 

	fi
}

condaCheck