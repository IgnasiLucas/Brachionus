#!/bin/bash

#Assemble expressed genes and transcripts

#B.plicatilis.fa is the reference genome
#B.plicatilis is the prefix for index files
#maker2.all.gff is the structural annotation

REFDIR=/data/eva/Brachionus/data/Transcriptoma/Genome
DATA=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation
TOPHAT=/home/eva/Brachionus/results/2016-12-19
CUFFLINKS=/home/eva/Brachionus/results/2017-01-24
CUFFMERGE=/home/eva/Brachionus/results/2017-01-24/cuffmerge


#Cuffdiff-quantifies expresSion

  cuffdiff -p 24 \
           -b $REFDIR/Brachionus_plicatilis_scaffold_min500.fasta \
      	   -u \
	   -L I1,P1,I2,P2 \
	   -q \
	   -o diff_out $CUFFMERGE/merged.gtf \
	   $TOPHAT/1A_S8/accepted_hits.bam,$TOPHAT/2A_S7/accepted_hits.bam,$TOPHAT/4A_S6/accepted_hits.bam \
	   $TOPHAT/3A_S9/accepted_hits.bam,$TOPHAT/5A_S2/accepted_hits.bam,$TOPHAT/6A_S3/accepted_hits.bam \
	   $TOPHAT/1C_S1/accepted_hits.bam,$TOPHAT/2C_S5/accepted_hits.bam,$TOPHAT/4C_S12/accepted_hits.bam \
	   $TOPHAT/3C_S11/accepted_hits.bam,$TOPHAT/5C_S4/accepted_hits.bam,$TOPHAT/6C_S10/accepted_hits.bam \
		1> cuffdiff.log 2> cuffdiff.error

