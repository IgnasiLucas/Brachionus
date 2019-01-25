#!/bin/bash

#B.plicatilis.fa is the reference genome
#B.plicatilis is the prefix for index files
#maker2.all.gff is the structural annotation

REFDIR=/data/eva/Brachionus/data/Transcriptoma/Genome
DATA=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation
TOPHAT=/home/eva/Brachionus/results/2016-12-19

#Cufflink-Assemble expressed genes and transcripts

for i in $DATA/*.fastq.gz; do
   echo $i
   dirname=$(basename "$i" | cut -d "." -f 1) #sería el mismo que name
   echo $dirname
   if [ ! -d $dirname.cufflinks ]; then
      mkdir $dirname.cufflinks #output assemble reads directory
   fi
	if [ ! -e $dirname.cufflinks/transcripts.gtf ]; then 
       cufflinks -p 24 \
                 -v \
                 -N \
		 -I 15000 \
		 --min-intron-length 10 \
                 -G $ANNOTATION/maker2.all.gff \
                 -o ./$dirname.cufflinks $TOPHAT/$dirname/accepted_hits.bam 1> $dirname.cufflinks/log 2> $dirname.cufflinks/error
   	fi
done

# Create a file called assemblies.txt that list assembly file for each sample

	if [ ! -e assemblies.txt ]; then
           find . -path "*/transcripts.gtf" > assemblies.txt #Para que liste todos los archivos
        fi

# Cuffmerge-Create a simple merged transcriptome annotation #input seria */transcripts.gtf

   if [ ! -e merged.gff ]; then
	cuffmerge -p 24 \
		  -g $ANNOTATION/maker2.all.gff \
                  -s $REFDIR/Brachionus_plicatilis_scaffold_min500.fasta \
		  -o ./ ./assemblies.txt 1> cuffmerge.log 2> cuffmerge.error
   fi

#Cudiff-quantifies expresion
