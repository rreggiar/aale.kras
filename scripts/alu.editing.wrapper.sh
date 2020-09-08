#!/usr/bin/bash

# wrapper script to generate alu loci for editing

input.counts.file=$1
output.tmp.counts=$2
output.locus.bed=$3

echo "$input.counts.file" "$output.tmp.counts" | Rscript alu.loci.for.editing.Rscript

until [ -f "$output.tmp.counts" ]

do
	sleep 5
done

echo "creating bed file"

cut -d"'" -f1 "$output.tmp.counts" | cut -d'=' -f2 | sed 's/:/\t/g ; s/-/\t/g ; s/_5//g' | sort -k1,1 -2,2n > "$output.locus.bed"

rm "$output.tmp.counts"
