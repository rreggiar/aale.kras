#!/bin/bash

inputdir=$1

for subdir in $inputdir/*/star.out/ ; do
  cd $subdir
  #sample=$(echo $PWD | cut -d'/' -f 9)
  sample=$(basename $subdir)
  echo $sample
  #exit
  #cd rna.edit.star.out
  #cd star.out
  #cd pass.2
  samtools view -S -b Aligned.out.sam > $sample.out.bam
  samtools sort -m 4G $sample.out.bam -o $sample.sorted.bam
  samtools index $sample.sorted.bam
  #bamCoverage -b $sample.sorted.bam -o $sample.bw
  #bamCoverage --filterRNAstrand forward -p 16\
  #  -b $sample.sorted.bam -o fwd.$sample.bw
  #bamCoverage --filterRNAstrand reverse -p 16\
  #  -b $sample.sorted.bam -o rev.$sample.bw


done
