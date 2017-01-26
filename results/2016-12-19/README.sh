#!/bin/bash

#B.plicatilis.fa is the reference genome
#B.plicatilis is the prefix for index files
#maker2.all.gff is the structural annotation

REFDIR=/data/eva/Brachionus/data/Transcriptoma/Genome
DATA=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation

if [ ! -d index ]; then mkdir index; fi

if [ ! -e index/B.plicatilis.fa ]; then
   ln -s $REFDIR/Brachionus_plicatilis_scaffold_min500.fasta index/B.plicatilis.fa
fi

if [ ! -e index/B.plicatilis.1.bt2 ]; then
   bowtie2-build index/B.plicatilis.fa index/B.plicatilis
fi

for i in $DATA/*.fastq.gz; do #Path para que vaya a los reads
   echo $i
   name=$(basename "$i" | cut -d "." -f 1) #name of each sample
   if [ ! -d $name ]; then
      mkdir $name #output alingments directory
   fi
   if [ ! -e $name.fastq ]; then
      gunzip -c $DATA/$name.trim.fastq.gz > $name.fastq
   fi

   #comprobar que options no dejar por defecto (ej. -p, etc)
   #min-intron-length=10 (5 in positive strand and 15 in negative strand); max-intron-length=15000 (respect to negative strand)
   if [ ! -e $name/accepted_hits.bam ]; then
      tophat2 -N 1 \
              --read-gap-length 1 \
              --read-realign-edit-dist 0 \
              -p 24 \
              -i 10 \
              -I 15000 \
              --max-multihits 1 \
	      -G $ANNOTATION/maker2.all.gff \
              -o ./$name ./index/B.plicatilis $name.fastq 1> $name/log 2> $name/error
      if [ -e $name/accepted_hits.bam ]; then
         rm $name.fastq
      fi
   fi
done


