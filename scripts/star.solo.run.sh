#!/bin/bash

INPUTDIR='/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/input/'
OUTPUTDIR='/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/output/star.solo'
GENOMEDIR='/public/groups/kimlab/genomes.annotations/HG.38/'

STAR --version

cd $INPUTDIR

STAR --genomeDir $GENOMEDIR \
  --runThreadN 12 \
  --readFilesIn $INPUTDIR/kras.fastq/*R2*.fastq $INPUTDIR/kras.fastq/*R1*.fastq \
  --soloOutFileNames solo.out \
  --outSAMattributes NH HI nM AS CB GN sM \
  --outSAMtype BAM SortedByCoordinate \
  --soloBarcodeReadLength 101 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist $INPUTDIR/10x.whitelist/737K-august-2016.txt 
