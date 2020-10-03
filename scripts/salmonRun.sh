#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

## activate correct env
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

