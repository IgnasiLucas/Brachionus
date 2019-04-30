#!/bin/bash
#
#				2019-03-29
#				==========
#
# The differential gene expression analysis was performed originall by Eva on
# 2017-01-31 with cuffdiff, which is quite standard. I would like here to try
# a different program. I choose EdgeR because seeems sofisticate and flexible,
# and because it can be installed with conda. To use it, I need count data:
# number of reads overlapping the transcripts, rather than normalized levels
# of expression such as FPKM. To get counts, I want to use RSEM (Li & Dewey, 2011).
# Unfortunately, I faild to install RSEM from conda, so that it's not available
# in the spec-file. I installed version 1.3.1 from source: https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz

MERGED=../2019-03-21/merged.gtf
REFERENCE=../../data/reference.fa
# BAMDIR=../2019-03-19
FASTQDIR=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
SAMPLE=( 1A_S8 1C_S1  2A_S7 2C_S5 3A_S9 3C_S11
         4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )

# It turns out that there is one exon in merged.gtf that does not have strand information.
# This makes rsem-prepare-reference halt. The exon is expressed in only one sample, 6C.
# But I don't want to just remove that record from merged.gtf. In fact, the original
# 6C_S10/transcripts.gtf includes the strand information, and I don't know why cuffmerge
# substituted the '+' sign for a '.' there. I feel compeled to re-write merged.gtf
# after correcting this type of problems

if [ ! -e z1.gtf ]; then
   cp $MERGED z1.gtf
   for exon in $(gawk '($7=="."){print $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5 "\\t"}' z1.gtf); do
      STRAND=$(find ../2019-03-21 -name transcripts.gtf -exec grep -P $exon '{}' \; | \
      gawk '{F[$7]++}END{for (f in F) print f "\t" F[f]}' | \
      sort -nrk 2,2 | head -n 1 | cut -f 1)
      echo $exon
      echo $STRAND
      sed -i -r "s/($exon[^\t]+\t)\.(\t.+)/\1\\$STRAND\2/" z1.gtf
   done
fi

if [ ! -d rsem-out ]; then mkdir rsem-out; fi

if [ ! -e rsem-out/merged.grp ]; then
   rsem-prepare-reference --gtf z1.gtf \
                          --bowtie \
                          --num-threads 24 \
                          $REFERENCE \
                          rsem-out/merged
fi

# We already aligned the reads with Tophat in 2019-03-19. Thus, it makes sense to take advantage
# of the 'convert-sam-for-rsem' script and save the alignment step here.
#
# for i in ${SAMPLE[@]}; do
#   if [ ! -d $i ]; then mkdir $i; fi
#   if [ ! -e $i/$i.bam ]; then
#      convert-sam-for-rsem --num-threads 4 \
#                           --memory-per-thread 2G \
#                           $BAMDIR/$i/accepted_hits.bam \
#                           $i/$i &
#   fi
# done
# wait
#
# Unfortunately, the original mapping was to the reference genome, not to the set of transcripts.
# Thus, it is not valid, and I need to let RSEM do the mapping.

for i in ${SAMPLE[@]}; do
   if [ ! -d $i ]; then mkdir $i; fi
   if [ ! -e $i/$i.isoforms.results ]; then
      rsem-calculate-expression --strandedness reverse \
                                --num-threads 4 \
                                --seed 19770319 \
                                --calc-pme \
                                --calc-ci \
                                --no-bam-output \
                                --fragment-length-mean 310.0 \
                                --fragment-length-sd 80.0 \
                                $FASTQDIR/$i.trim.fastq.gz \
                                rsem-out/merged \
                                $i/$i &
   fi
   if [ ! -e $i/$i.pdf ]; then
      rsem-plot-model $i/$i $i/$i.pdf
   fi
done
wait

# Below I extract the counts from every sample and paste them together in one file
# for genes and one for isoforms. I checked that results from all samples have the
# same number of rows, and in the same order:
#
# find . -name '*.genes.results' -exec bash -c 'cut -f 1 $1 | md5sum' _ '{}' \;
# find . -name '*.isoforms.results' -exec bash -c 'cut -f 1 $1 | md5sum' _ '{}' \;

COUNTS=( zero one two three four ExpCount six seven PostCount PostCount)
WHAT=( genes isoforms )
HOW=( ExpCount PostCount )
COLUMN=( 5 5 8 9)
#        0 1 2 3  <-  echo "$what * 2^0 + $how * 2^1" | bc -l
# what   0 1 0 1
# how    0 0 1 1
#
# In RSEM results files, the expected counts are in the fifth column, both for genes
# and for isoforms. But, the posterior mean counts are in the eigth column for genes
# and in the nineth, for isoforms.
for what in 0 1; do
   for how in 0 1; do
      if [ ! -e ${WHAT[$what]}.${HOW[$how]}.txt ]; then
         echo -e "${WHAT[$what]}\t1A_S8\t1C_S1\t2A_S7\t2C_S5\t3A_S9\t3C_S11\t4A_S6\t4C_S12\t5A_S2\t5C_S4\t6A_S3\t6C_S10" > ${WHAT[$what]}.${HOW[$how]}.txt
         INDEX=$( echo "$what * 2^0 + $how * 2^1" | bc -l )
         paste <(cut -f 1,${COLUMN[$INDEX]} 1A_S8/1A_S8.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   1C_S1/1C_S1.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   2A_S7/2A_S7.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   2C_S5/2C_S5.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   3A_S9/3A_S9.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]} 3C_S11/3C_S11.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   4A_S6/4A_S6.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]} 4C_S12/4C_S12.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   5A_S2/5A_S2.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   5C_S4/5C_S4.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]}   6A_S3/6A_S3.${WHAT[$what]}.results | tail -n +2) \
               <(cut -f ${COLUMN[$INDEX]} 6C_S10/6C_S10.${WHAT[$what]}.results | tail -n +2) \
               >> ${WHAT[$what]}.${HOW[$how]}.txt
      fi
   done
done

if [ ! -e Expected_posterior.png ]; then
   Rscript plot.R
fi

# CONCLUSION
# ----------
#
# Both the expected and the posterior mean count are equivalent.
