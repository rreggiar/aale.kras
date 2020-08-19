#!/bin/bash

INPUTDIR='/public/home/rreggiar/projects/aale.kras/data/single.cell.rna.seq/input'
OUTPUTDIR='/public/home/rreggiar/projects/aale.kras/data/single.cell.rna.seq/output'
GCINDEX='/public/groups/kimlab/indexes/gencode.32.v.1.index'
TEINDEX='/public/groups/kimlab/indexes/gencode.te.locus.v.1.index/'
TX2GENE='/public/groups/kimlab/genomes.annotations/gen.32.ucsc.rmsk.tx.2.gene.tsv'

salmon -v

cd $INPUTDIR

for INDEX in $TEINDEX ; do

		salmon alevin -i $INDEX \
			-l ISR \
			-1 ${INPUTDIR}/'kras.fastq/Kim_AL_K_Stable_S2_L001_R1_001-002.fastq'\
			-2 ${INPUTDIR}/'kras.fastq/Kim_AL_K_Stable_S2_L001_R2_001-003.fastq'\
      --chromium \
			-p 16 \
			--tgMap $TX2GENE \
      --o $OUTPUTDIR/$(basename $INDEX).alevin.out
    
done

