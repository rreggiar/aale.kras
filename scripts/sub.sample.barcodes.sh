#!/bin/bash

INPUTDIR='/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/solo.outGene/'
OUTPUTDIR='/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/cluster.bam.files/'
# establish the constant split goal
num_files=3

cd $INPUTDIR
mkdir splits.out
# iterate over the barcode files I've exported from seurat
for BCFILE in *.barcodes.txt; do
  # number of cells per file (by barcode)
  num_cells=$(wc -l <$BCFILE)
  # calculate how many cells to include in each split
  ((cells_per_file=(num_cells + num_files - 1) / num_files))
  # split into 3 files with the leader $BCFILE.split.
  split --lines=$cells_per_file $BCFILE $BCFILE.split.

  # iterate over splits and use them to subset the bam
  for SUB in $BCFILE.split.*; do
    samtools view -O BAM -@ 12 \
      -D "CB:$SUB" \
      -o $OUTPUTDIR/$SUB.bam \
      Aligned.sortedByCoord.out.bam
    # dump the splits in a single file
    mv $SUB splits.out
  done
done
