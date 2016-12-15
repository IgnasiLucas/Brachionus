#!/bin/bash
#
#				2016-12-14
#				----------
#
# I want to know the distributions of lengths and qualities of the reads.
# This information seems to be available for the original reads, in the
# zipped fastqc reports. But I want to look at the trimmed reads.

FASTQDIR=/data/eva/Brachionus/data/Transcriptoma/Trim_reads
SAMPLE=(
   1A_S8
   1C_S1
   2A_S7
   2C_S5
   3A_S9
   3C_S11
   4A_S6
   4C_S12
   5A_S2
   5C_S4
   6A_S3
   6C_S10)

for i in "${SAMPLE[@]}"; do
   if [ ! -e $i.stats ]; then
      calculate_stats -o $i.stats $FASTQDIR/$i.trim.fastq.gz &
   fi
done
wait

if [ ! -e summary_fastq.txt ]; then
   echo -e "Sample\tMinLength\tAverage\tMaxLength\tNumSeqs\tQ20\tQ30\tMinQ\tAverage\tMaxQ" > summary_fastq.txt
   for i in "${SAMPLE[@]}"; do
      gawk -v SAMPLE=$i 'BEGIN{
         LENGTH=1
      }(/^Quality stats and distribution/){
         LENGTH=0
      }((/^minimum:/) && (LENGTH == 1)){
         MINLEN = $2
      }((/^maximum:/) && (LENGTH == 1)){
         MAXLEN = $2
      }((/^average:/) && (LENGTH == 1)){
         AVELEN = $2
      }(/^num. seqs.:/){
         NUMSEQ = $3
      }(/^Q20:/){
         Q20 = $2
      }(/^Q30:/){
         Q30 = $2
      }((/^minimum:/) && (LENGTH == 0)){
         MINQ = $2
      }((/^maximum:/) && (LENGTH == 0)){
         MAXQ = $2
      }((/^average:/) && (LENGTH == 0)){
         AVEQ = $2
      }END{
         print SAMPLE "\t" MINLEN "\t" AVELEN "\t" MAXLEN "\t" NUMSEQ "\t" Q20 "\t" Q30 "\t" MINQ "\t" AVEQ "\t" MAXQ
      }' $i.stats >> summary_fastq.txt
   done
fi

# CONCLUSIONS
# -----------
#
# Trimmed reads are between 20 and 75 bp, and the average is >74 bp in all samples.
# Average qualities ~34.
#
# I used calculate_stats from seqcrumbs, which is not as fast as vsearch --fastq_stats.
# However, the output is easier to parse.
