#!/bin/bash
#
#				2019-03-19
#				==========
#
# Here I want to replicate the mapping of reads with tophat. I notice that originally
# we did not specify the library type. Note that I use a different conda environment
# to run this analysis, which requires tophat. The environment is saved in the local
# spec-file.txt.

REFERENCE=/data/eva/Brachionus/data/Transcriptoma/Genome/Brachionus_plicatilis_scaffold_min500.fasta
FASTQDIR=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation/maker2.all.gff
SAMPLE=( 1A_S8 1C_S1  2A_S7 2C_S5 3A_S9 3C_S11
         4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )

if [ ! -d index ]; then mkdir index; fi

if [ ! -e index/B.plicatilis.fa ]; then
   ln -s $REFERENCE index/B.plicatilis.fa
fi

if [ ! -e index/B.plicatilis.1.bt2 ]; then
   bowtie2-build index/B.plicatilis.fa index/B.plicatilis
fi

for name in ${SAMPLE[@]}; do
   if [ ! -d $name ]; then
      mkdir $name #output alingments directory
   fi

   if [ ! -e $name/accepted_hits.bam ]; then
      tophat  --read-mismatches 1 \
              --read-gap-length 1 \
              --read-realign-edit-dist 0 \
              --min-intron-length 10 \
              --max-intron-length 15000 \
              --max-multihits 1 \
              --library-type fr-firststrand \
              --num-threads 24 \
	      --GTF $ANNOTATION \
              --rg-id $name \
              --rg-sample $name \
              --output-dir ./$name \
              ./index/B.plicatilis $FASTQDIR/$name.trim.fastq.gz 1> $name/log 2> $name/error
   fi
done

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
      cut -f 5 $i/deletions.bed | tail -n +2 | sed 's/^/-/' > zindels
      cut -f 5 $i/insertions.bed | tail -n +2 >> zindels
      R --no-save --args $i < plot.R
      mv indels.png $i/indels.png
      rm zindels
   fi
done

# The junctions.bed file describes the coordinates of the reference genome that
# are connected by adjacent exons. That is, the introns. The very last number of
# each line is the length of the intron. It would be nice to make sure that the
# distribution of intron lengths are the same across samples. This time, I make
# only one plot with all samples in it.

if [ ! -e introns.png ] || [ $(find . -name junctions.bed -newer introns.png) -gt 0 ]; then
   for i in ${SAMPLE[@]}; do
      tail -n +2 $i/junctions.bed | cut -f 12 | cut -d ',' -f 2 > z$i
   done
   R --no-save < plot2.R
   rm z*
fi

# CONCLUSIONS
# ===========
#
# The small modification of the arguments passed to Tophat did not result in any
# difference in the mapping, as assessed by the mapping success and the distribution
# of indel lengths. All samples have equivalent intron length distributions, which
# confirms the homogeneous quality of the samples.

# CLEAN UP
# ========

for i in ${SAMPLE[@]}; do
   if [ -d $i/logs ];              then rm -r $i/logs; fi
   if [ -e $i/align_summary.txt ]; then rm $i/align_summary.txt; fi
   if [ -e $i/deletions.bed ];     then rm $i/deletions.bed; fi
   if [ -e $i/error ];             then rm $i/error; fi
   if [ -e $i/insertions.bed ];    then rm $i/insertions.bed; fi
   if [ -e $i/junctions.bed ];     then rm $i/junctions.bed; fi
   if [ -e $i/prep_reads.info ];   then rm $i/prep_reads.info; fi
   if [ -e $i/unmapped.bam ];      then rm $i/unmapped.bam; fi
done
