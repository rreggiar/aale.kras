#!/bin/bash

INPUTDIR='/public/groups/kimlab/aale.kras/data/bulk.rna.seq/exo/input'
#OUTDIR=quantFiles
TXINDEX='/public/groups/kimlab/indexes/gencode.32.v.1.index/'
#TXINDEX='/public/groups/kimlab/indexes/te.locus.v.1.index/'
HG38='/public/groups/kimlab/genomes.annotations/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
GENOMEDIR='/public/groups/kimlab/genomes.annotations/HG.38.w.te.cons/'
ADAPTERS='/public/groups/kimlab/genomes.annotations/adapters'
OUTDIR=gencode.salmon.out

STAR --version
salmon -v
trimmomatic -version

cd $INPUTDIR
#mkdir $OUTDIR

for SAMPLE in $PWD/* ; do
    echo $SAMPLE

	cd $SAMPLE
	NAME=$(basename $SAMPLE)

	if [ ! -f *_paired.fq.gz ]; then
	
		read1=*_R1_001.fastq.gz
		read2=*_R2_001.fastq.gz
	    echo "Read 1:" ${read1}
		echo "Read 2:" ${read2}
        echo "fastq files must be trimmed"
        echo 'Trimming Now..'
        trimmomatic PE -threads 12 ${read1} ${read2} output_forward_paired.fq.gz \
            output_forward_unpaired.fq.gz output_reverse_paired.fq.gz \
            output_reverse_unpaired.fq.gz \
            ILLUMINACLIP:$ADAPTERS/NexteraPE-PE.fa:1:30:10:4:true LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36
    fi

	if [ -f $SAMPLE/$OUTDIR/quant.sf ]; then
		echo "SALMON HAS ALREADY BEEN RUN"

	else
		echo "SALMON HAS NOT BEEN RUN, NOW RUNNING" ${read1}
		
		mkdir $SAMPLE/$OUTDIR/
		
        trim1=$SAMPLE/*output_forward_paired.fq.gz
        trim2=$SAMPLE/*output_reverse_paired.fq.gz

		salmon quant -i ${TXINDEX} \
			--libType A \
			-1 ${trim1}\
			-2 ${trim2}\
			-p 64\
			--validateMappings \
			--gcBias \
			--seqBias \
			--output $SAMPLE/$OUTDIR/
    fi
    '''
    echo "SALMON COMPLETE, CHECKING STAR"
    trim1=$SAMPLE/output_forward_paired.fq.gz
    trim2=$SAMPLE/output_reverse_paired.fq.gz
    if [ ! -f "$SAMPLE/rna.edit.star.out/Aligned.out.sam" ]; then
        echo "RUNNING 1 PASS STAR for editing ON" $SAMPLE
        firstRunDir=$SAMPLE/rna.edit.star.out/
        mkdir $firstRunDir
        cd $firstRunDir
        #secondRunDir=$SAMPLE/star.out/pass.2/
        #mkdir $secondRunDir
        STAR --genomeDir $GENOMEDIR \
        --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) \
        --outFilterMatchNminOverLread 0.95 \
        --outSAMtype BAM SortedByCoordinate
    fi   
    echo "SALMON COMPLETE, CHECKING STAR"
    trim1=$SAMPLE/output_forward_paired.fq.gz
    trim2=$SAMPLE/output_reverse_paired.fq.gz
    if [ ! -f *Aligned.out.sam ]; then
        echo "RUNNING 2 PASS STAR ON" $SAMPLE
        firstRunDir=$SAMPLE/rna.edit.star.out/
        mkdir $firstRunDir
        cd $firstRunDir
        #secondRunDir=$SAMPLE/star.out/pass.2/
        #mkdir $secondRunDir
        STAR --genomeDir $GENOMEDIR --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) --outFilterMatchNminOverLread 0.95
    fi
done
cd $INPUTDIR
for SAMPLE in $PWD/*; do
    trim1=$SAMPLE/output_forward_paired.fq.gz
    trim2=$SAMPLE/output_reverse_paired.fq.gz
    cd $SAMPLE/star.out/pass.2/
	STAR --genomeDir $GENOMEDIR --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) --sjdbFileChrStartEnd $INPUTDIR/*/star.out/SJ.out.tab
    '''
done

