#!/bin/bash

#B.plicatilis.fa is the reference genome
#B.plicatilis is the prefix for index files


REFDIR=/data/eva/Brachionus/data/Transcriptoma/Genome
DATA=/data/eva/Brachionus/data/Transcriptoma/Trim_reads

if [ ! -d index ]; then mkdir index; fi

if [ ! -e index/B.plicatilis.fa ]; then
   ln -s $REFDIR/B_plicat_PLATANUS.fasta index/B.plicatilis.fa
fi

if [ ! -e index/B.plicatilis.1.bt2 ]; then
   bowtie2-build index/B.plicatilis.fa index/B.plicatilis #Falta poner path (the recommendatio:all files in he same directory)
fi

for i in $DATA/*.fastq.gz; do #Path para que vaya a los reads
   echo $i
   name = $(basename "$i" | cut -d "." -f 1) #name of each sample
   if [ ! -d $name ]; then
      mkdir $name #output alingments directory
   fi
   if [ ! -e $name.fastq ]; then
      gunzip -c $DATA/$name_.fastq.gz > $name.fastq
   fi

   #comprobar que options no dejar por defecto (ej. -p, etc)
   if [ -e $name/accepted_hits.bam ]; then
      tophat2 -N 1 \
              --read-gap-length 1 \
              --read-realign-edit-dist \
              -p 24 \
              -i    \
              -I    \
              --max-multihits 1 \
              -o ./$name ./index/B.plicatilis $name.fastq 1> $name/log 2> $name/err

      rm $name.fastq
   fi
done


