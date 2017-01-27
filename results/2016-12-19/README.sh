#!/bin/bash

#B.plicatilis.fa is the reference genome
#B.plicatilis is the prefix for index files
#maker2.all.gff is the structural annotation

REFDIR=/data/eva/Brachionus/data/Transcriptoma/Genome
DATA=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation
SAMPLE=( 1A_S8 1C_S1  2A_S7 2C_S5 3A_S9 3C_S11
         4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )

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

   if [ ! -e $name/accepted_hits.bam ]; then
      if [ ! -e $name.fastq ]; then
         gunzip -c $DATA/$name.trim.fastq.gz > $name.fastq
      fi

   #comprobar que options no dejar por defecto (ej. -p, etc)
   #min-intron-length=10 (5 in positive strand and 15 in negative strand); max-intron-length=15000 (respect to negative strand)
      tophat2 -N 1 \
              --read-gap-length 1 \
              --read-realign-edit-dist 0 \
              -p 24 \
              -i 10 \
              -I 15000 \
              --max-multihits 1 \
	      -G $ANNOTATION/maker2.all.gff \
              -o ./$name ./index/B.plicatilis $name.fastq 1> $name/log 2> $name/error

      if [ -e $name/accepted_hits.bam ] && [ -e $name.fastq ]; then
         rm $name.fastq
      fi
   fi
done

# I want to create the summary table if it doesn't exist or if it
# is older than any algin_summary.txt file in this folder.
if [ ! -e summary_mapping.txt ] || [ $(find . -name align_summary.txt -newer summary_mapping.txt | wc -l) -gt 0 ]; then
   echo -e "Sample\tReads\tMapped\tRate" > summary_mapping.txt
   for i in ${SAMPLE[@]}; do
      gawk -v SAMPLE=$i '(/Input/){
         READS = $3
      }(/Mapped/){
         MAPPED = $3
      }END{
         RATE = 100 * MAPPED / READS
         printf("%s\t%u\t%u\t%.2f \%\n", SAMPLE, READS, MAPPED, RATE)
      }' $i/align_summary.txt >> summary_mapping.txt
   done
fi

# The insertions.bed and deletions.bed files of each sample have a lot of information
# that may be worth summarizing with a simple histogram of sizes of insertions and
# deletions. We could make a single figure with all samples, but to keep the script
# simple, I will do one figure for each sample.

for i in ${SAMPLE[@]}; do
   if [ ! -e $i/indels.png ]; then
      # First I create a temporal file with the data of this sample
      # to pass to R. It will be just one list of sizes, where deletions
      # have negative size:
      cut -f 5 $i/deletions.bed | tail -n +2 | sed 's/^/-/' > zindels
      cut -f 5 $i/insertions.bed | tail -n +2 >> zindels
      # This makes a plot in this folder, and I move it to its sample
      # folder. Even though I'm using --args to pass the name of the sample
      # to the R script, I find it easier to move the figure after it's done
      # than to tell the script how to use the name to save it in the right folder.
      R --no-save --args $i < plot.R
      mv indels.png $i/indels.png
      rm zindels
   fi
done

# The junctions.bed file describes the coordinates of the reference genome that
# are connected by adjacent exons. That is, the introns. The very last number of
# each line is the length of the intron. It would be nice to make sure that the
# distribution of intron lengths are the same across samples. This time, I make
# only one plot with all samples in it. Actually, I don't have time to do this...

#if [ ! -e introns.png ] || [ $(find . -name junctions.bed -newer introns.png) -gt 0 ]; then
#   for i in ${SAMPLE[@]}; do
#      tail -n +2 $i/junctions.bed | cut -f 12 | cut -d ',' -f 2 > z$i
#   done
#   R --no-save < plot2.R
#fi
