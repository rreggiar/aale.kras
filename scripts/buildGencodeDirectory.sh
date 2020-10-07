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
gencodeTransctiptFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".transcripts.fa.gz"
gencodePrimaryAssemblyFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".transcripts.fa.gz"
gencodeLncRNATranscriptFA="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"$version"/gencode.v"$version".lncRNA_transcripts.fa.gz"
ucscRmskInsertFA="/public/groups/kimlab/genomes.annotations/formatted.UCSC.gb.rmsk.insert.fa"

# generate destination directory
if [ ! -d "$outputDir" ]; then

	set -x
	mkdir "$outputDir"
	set +x

fi

# symlink in the rmsk reference -- ln won't overwrite by default
ln -s "$ucscRmskInsertFA" "$outputDir"

# for succint iteration
dataGenerationList=("$gencodeAnnotationGTF" "$gencodeTransctiptFA" "$gencodePrimaryAssemblyFA" "$gencodeLncRNATranscriptFA")

function downloadDataSets(){

	dataGenerationList="$1"

	# non-redundant download of necc. data for annotation directory construction
	for targetData in "$dataGenerationList"; do

		destinationPath="$outputDir"/"$(basename "$targetData")"

		set -x
		if [ ! -f "$destinationPath" ]; then

			wget "$targetData" --output-document="$destinationPath"

		fi
		set +x

	done

}


downloadDataSets "$dataGenerationList"











