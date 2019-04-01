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
BAMDIR=../2019-03-19
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
fi

for exon in $(gawk '($7=="."){print $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5 "\\t"}' z1.gtf); do
   STRAND=$(find ../2019-03-21 -name transcripts.gtf -exec grep -P $exon '{}' \; | \
   gawk '{F[$7]++}END{for (f in F) print f "\t" F[f]}' | \
   sort -nrk 2,2 | head -n 1 | cut -f 1)
   echo $exon
   echo $STRAND
   sed -i -r "s/($exon[^\t]+\t)\.(\t.+)/\1\\$STRAND\2/" z1.gtf
done

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

for i in ${SAMPLE[@]}; do
   if [ ! -d $i ]; then mkdir $i; fi
   if [ ! -d $i/$i.bam ]; then
      convert-sam-for-rsem --num-threads 4 \
                           --memory-per-thread 2G \
                           $BAMDIR/$i/accepted_hits.bam \
                           $i/$i &
   fi
done
wait

for i in ${SAMPLE[@]}; do
   if [ ! -e $i/$i.isoforms.results ]; then
      rsem-calculate-expression --alignments \
                                --strandedness reverse \
                                --num-threads 4 \
                                --seed 19770319 \
                                --calc-pme \
                                --calc-ci \
                                --no-bam-output \
                                --fragment-length-mean 200.0 \
                                --fragment-length-sd 80.0 \
                                $i/$i.bam \
                                rsem-out/merged
   fi
done
wait

