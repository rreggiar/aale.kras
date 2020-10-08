#!/bin/bash

# rreggiar@ucsc.edu

# download gencode materials for a provided release and construct materials used in standard analytical procedures:
## tx2gene -- transcript to gene annotation that is used by DESeq
## pureLncRna -- genes that have only non-coding isoforms
## ucsc.rmsk.* -- combining the gencode annotation with the repeat masked 'transcripts' from UCSC GB

# generate a couple of salmon indexes that will enable standard kimlab quantification procedures
## [selective alignment](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
## [process-aware alignmnet (modified below)](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

# gencode data links
## [gencode data summary page](https://www.gencodegenes.org/human/)
## [gencode comprehensive annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".annotation.gtf.gz)
## [gencode transcript sequences](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".transcripts.fa.gz)
## [gencode primary assembly](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/GRCh38.primary_assembly.genome.fa.gz)
## [gencode lncRNA transcript sequences](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".lncRNA_transcripts.fa.gz)

# ucsc genome browser
## [ucsc rmsk insertion fasta](/public/groups/kimlab/genomes.annotations/formatted.UCSC.gb.rmsk.insert.fa)
### replace all the ' ' with '_' so that it plays nicely with downstream processes
## [ucsc repeat browser consensus fasta](/public/groups/kimlab/genomes.annotations/repeat.browser.hg38reps.fa)
### NOTE: these used to be accessible via the UCSC GB FTP server but I cannot locat them anymore ... as always they can be 
### generated via the `table browser` utility but this does not support programmatic access

scriptName=$(basename $0)
if [ $# -lt 1 ]; then
    echo "error: usage $scriptName  gencode_version"
    echo "example: $scriptName 35"
    exit 1
fi

command="$@"
echo "$scriptName" "$command"

# references to enable download and placement of datasets
kimlabGenomesDir="/public/groups/kimlab/genomes.annotations"
kimlabIndexDir="/public/groups/kimlab/indexes"
version="$1"
outputDir="$kimlabGenomesDir"/"gencode.""$version"

# gencode-format for FTP download, they keep everything nice and regular just modifying the version number
gencodeAnnotationGTF="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".annotation.gtf.gz"
gencodeTranscriptFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".transcripts.fa.gz"
gencodePrimaryAssemblyFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/GRCh38.primary_assembly.genome.fa.gz"
gencodeLncRNATranscriptFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".lncRNA_transcripts.fa.gz"
ucscRmskInsertFA="/public/groups/kimlab/genomes.annotations/formatted.UCSC.gb.rmsk.insert.fa"
ucscRmskInsertTx2GeneCSV="/public/groups/kimlab/genomes.annotations/formatted.UCSC.gb.rmsk.insert.tx.to.gene.csv"

# generate destination directory
if [ ! -d "$outputDir" ]; then

	set -x
	mkdir "$outputDir"
	set +x

fi

# symlink in the rmsk reference -- ln won't overwrite by default
ln -s "$ucscRmskInsertFA" "$outputDir"
# and the rmsk tx2gene
ln -s "$ucscRmskInsertTx2GeneCSV" "$outputDir"

# for succint iteration
dataGenerationList=("$gencodeAnnotationGTF" "$gencodeTranscriptFA" "$gencodePrimaryAssemblyFA" "$gencodeLncRNATranscriptFA")

function downloadDataSets(){

	dataGenerationList="$1"

	# non-redundant download of necc. data for annotation directory construction
	for targetData in "${dataGenerationList[@]}"; do

		destinationPath="$outputDir"/"$(basename "$targetData")"

		set -x
		if [ ! -f "$destinationPath" ]; then

			wget "$targetData" --output-document="$destinationPath"

		fi
		set +x

	done

}

function makeTx2Gene(){

	gencodeTranscriptFA="$1"
	ucscRmskInsertTx2GeneCSV="$2"
	outputDir="$3"

	set -x

	# gencode tx names
	if [ ! -f "$outputDir"/"gencode.v"$version".transcript.names.txt" ]; then
		zcat "$outputDir"/"$(basename "$gencodeTranscriptFA")" | grep '>' | cut -d'>' -f2 > "$outputDir"/"gencode.v"$version".transcript.names.txt"
	fi	
	# gene names
	if [ ! -f "$outputDir"/"gencode.v"$version".gene.names.txt" ]; then
		cut -d'|' -f6 "$outputDir"/"gencode.v"$version".transcript.names.txt" > "$outputDir"/"gencode.v"$version".gene.names.txt"
	fi
	# tx to gene
	if [ ! -f "$outputDir"/"gencode.v"$version".tx.to.gene.csv" ]; then
		paste -d, "$outputDir"/"gencode.v"$version".transcript.names.txt" "$outputDir"/"gencode.v"$version".gene.names.txt" > "$outputDir"/"gencode.v"$version".tx.to.gene.csv"
	fi
	# gencode + rmsk tx.2.gene
	if [ ! -f "$outputDir"/"gencode.v"$version".ucsc.rmsk.tx.to.gene.csv" ]; then
		cat "$outputDir"/"gencode.v"$version".tx.to.gene.csv" "$ucscRmskInsertTx2GeneCSV" > "$outputDir"/"gencode.v"$version".ucsc.rmsk.tx.to.gene.csv"
	fi

	set +x

}

function makeSalmonDecoys(){

	version="$1"
	outputDir="$2"
	gencodePrimaryAssemblyFA="$3"
	decoysOut="$outputDir"/"gencode.v""$version"".decoys.txt"

		if [ ! -f "$decoysOut" ]; then

			set -x

			grep "^>" <(gunzip -c "$outputDir"/"$(basename "$gencodePrimaryAssemblyFA")") | \
				cut -d" " -f1 | \
				sed -e 's/>//g' > "$decoysOut"

			set +x
		fi

}

function makeProcessAwareReferences(){

	outputDir="$1"

	set -x

	echo "$outputDir" | Rscript processAwareSalmonReference.R

	echo "$outputDir" | Rscript genereLinkedTxome.R

	set +x

}

function makeSalmonIndexes(){

	if [ ! -x "$(command -v salmon)" ]; then

		echo "please install salmon or activate the correct conda environment"
		exit 1

	else

		set -x

		version="$1"
		outputDir="$2"
		kimlabIndexDir="$3"
		gencodeTranscriptFA="$4"
		ucscRmskInsertFA="$5"
		gencodePrimaryAssemblyFA="$6"
		processAwareFA="$outputDir"/"$7"


		genTxFA="$outputDir"/"$(basename "$gencodeTranscriptFA")"
		genomeFA="$outputDir"/"$(basename "$gencodePrimaryAssemblyFA")"

		salmonVersion=`salmon -v | cut -d' ' -f2`

		salmon index \
			-t <(cat "$genTxFA" "$genomeFA") \
			-i "$kimlabIndexDir"/"sel.align.gencode.v""$version"".salmon.v""$salmonVersion"".sidx" \
			-p 16 \
			-d "$outputDir"/"gencode.v""$version"".decoys.txt"

		
		salmon index \
			-t <(cat "$genTxFA" "$ucscRmskInsertFA" "$genomeFA") \
			-i "$kimlabIndexDir"/"sel.align.gencode.v""$version"".ucsc.rmsk.salmon.v""$salmonVersion"".sidx" \
			-p 16 \
			-d "$outputDir"/"gencode.v""$version"".decoys.txt"

		salmon index \
			-t <(cat "$processAwareFA" "$genomeFA") \
			-i "$kimlabIndexDir"/"sel.align.gencode.v""$version"".process.aware.salmon.v""$salmonVersion"".sidx" \
			-p 16 \
			-d "$outputDir"/"gencode.v""$version"".decoys.txt"

		set +x

	fi

}

downloadDataSets "$dataGenerationList"

makeTx2Gene "$gencodeTranscriptFA" "$ucscRmskInsertTx2GeneCSV" "$outputDir"

makeSalmonDecoys "$version" "$outputDir" "$gencodePrimaryAssemblyFA"

makeProcessAwareReferences "$outputDir"

makeSalmonIndexes "$version" "$outputDir" "$kimlabIndexDir" \
	"$gencodeTranscriptFA" "$ucscRmskInsertFA" "$gencodePrimaryAssemblyFA" \
	"gencode.v""$version"".annotation.expanded.fa"
