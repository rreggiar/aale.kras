#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env
# run trimmomatic to remove adapter sequences from the ends of reads
# the resulting paired output files will be used for all downstream analysis

scriptName=$(basename $0)
if [ $# -lt 2 ]; then
    echo "error: usage $scriptName  sample.input.dir adapterChoice"
    echo "example: $scriptName /path/to/{ctrl,kras}.{1,2,3..} {exo, intra}"
    exit 1
    # rerwip: Make the outDir argument make sense for symlinked databases
fi

# this is a constant 
adapterPath='/public/groups/kimlab/genomes.annotations/adapters'
# which adapters we use will depend on the library preparation
# the main situations this is relevant to are the different intracellular
# and extracellular preparations. This allows us to specify either of those as
# a shortcut or provide a custom argument that will be checked against our
# adapter library
if [[ "$2" == 'intra' ]]; then
	adapterChoice='TruSeq3-PE.fa'
elif [[ "$2" == 'exo' ]]; then
	adapterChoice='NexteraPE-PE.fa'
else
	if [[ -f "$adapterPath"/"$2" ]]; then
			adapterChoice="$2"
	else
		echo "you have provided an adapter file that doesn't exist (yet?): "$2""
		exit 1
	fi
fi

## don't want dateStamp functionality here
# dateStamp="$3"
# set -x
# echo "script: $scriptName"
# echo "time: $dateStamp"
# set +x 

inputDir="$1"

set -x
echo "input: $inputDir"
echo "trimmomatic: $(trimmomatic -version)"
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

## since trimming is deterministic and reproducible, I'm going to treat these as un-dated
## output files for now so they can be accessed by various tools/pipelines and not be needlessly
## reproduced

## currently, this runs assuming you want output sent to $PWD, as in we are iterating over
## input dir's. This is somthing I need to talk to aedavids about...
function runTrimmomatic() {
	inputDir="$1"

	# check for existence of canonical output file name
	if [ ! -f "$1"/output_forward_paired.fq.gz ]; then
		# if the read files exist and are symlinks
		if [[ -L "$1"/*_R1_001.fastq.gz ]]; then
			# assign the reads  -- prepping for move to $SCRATCH data directories
			read_1=$(readlink *_R1_001.fastq.gz)
			read_2=$(readlink *_R2_001.fastq.gz)
		else
			read_1=*_R1_001.fastq.gz
			read_2=*_R2_001.fastq.gz
		fi
		
		set -x
		echo "read 1:" "$read_1"
		echo "read 2:" "$read_2"
		echo "Trimming fastq"

		# run trimm with recommended arguments
		trimmomatic PE  -threads 8 \
		"$read_1" "$read_2" \
		output_forward_paired.fq.gz \
		output_forward_unpaired.fq.gz \
		output_reverse_paired.fq.gz \
		output_reverse_unpaired.fq.gz \
		ILLUMINACLIP:"$adapters"/"$lib_choice":1:30:10:4:true \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36
		set +x

	fi
}

condaCheck

runTrimmomatic "${inputDir}"

