#!/bin/bash
#
#				2019-03-21
#				==========
#
# Here I want to run cufflinks again, as Eva did in 2017-01-24, just to check
# the effect of some options that were not used back then, such as the library
# type. Also, the mean fragment length is 310, not 200 bp that is set by default.
# We do not know the standard deviation of the library size, so I leave it at 80.

REFERENCE=/data/eva/Brachionus/data/Transcriptoma/Genome/Brachionus_plicatilis_scaffold_min500.fasta
BAMDIR=../2019-03-19
ANNOTATION=/data/eva/Brachionus/data/Transcriptoma/Annotation/maker2.all.gff
SAMPLE=( 1A_S8 1C_S1  2A_S7 2C_S5 3A_S9 3C_S11
         4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )


#Cufflink-Assemble expressed genes and transcripts.
for i in ${SAMPLE[@]}; do
   if [ ! -d $i ]; then mkdir $i; fi
   if [ ! -e $i/transcripts.gtf ]; then
      cufflinks --output-dir $i \
                --num-threads 24 \
                --GTF-guide $ANNOTATION \
                --multi-read-correct \
                --library-type fr-firststrand \
                --frag-len-mean 310 \
                --max-intron-length 15000 \
                --min-intron-length 10 \
                $BAMDIR/$i/accepted_hits.bam 1> $i/cufflinks.log 2> $i/cufflinks.err
   fi
   # This will reduce the size of the cufflinks.err files from 13M to a few bytes:
   if [ -e $i/cufflinks.err ]; then
      gawk -v FS="\r" '{print $NF}' $i/cufflinks.err > z1
      mv z1 $i/cufflinks.err
   fi

   if [ ! -e $i/fpkm_hist.png ]; then
      gawk '{gsub(/;|\"/,""); for (i=9; i<=NF; i++) {if ($i == "FPKM") print $(i+1)}}' $i/transcripts.gtf > z1
      Rscript plot.R $i/fpkm_hist.png z1
      rm z1
   fi
   # After looking at the histograms of relative expression level, I notice that all samples have
   # a kind of long tail of a few highly expressed genes. I arbitrarily establish the threshold of
   # highly expressed genes as those with PFKM > 1000. Below, I extract a list of those transcripts
   # and exons for every sample, so that we can compare them among samples. If there was a common
   # set of highly expressed genes in all samples, I would consider removing them from analysis with
   # the --mask-file option to cufflinks.
   if [ ! -e $i/superexpressed.gtf ]; then
      gawk -v FS="\t" '{
         gsub(/;|\"/,"")
         split($9,A,/ /)
         for (i=1; i<=length(A); i++) {
            if ((A[i] == "FPKM") && (A[i+1] + 0 > 1000.0)) print $0
         }
      }' $i/transcripts.gtf > $i/superexpressed.gtf
   fi
done

# Compare the superexpressed genes among samples, and list the ones that are common to all them.
if [ ! -e common_superexpressed.gtf ]; then
   gawk '{F[$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8]++}END{
      for (f in F) {
         G[F[f]]++
         if (F[f] == 12) print f > "z1"
      }
      for (g in G) {
         print g "\t" G[g] > "z2"
      }
   }' 1A_S8/superexpressed.gtf \
      1C_S1/superexpressed.gtf \
      2A_S7/superexpressed.gtf \
      2C_S5/superexpressed.gtf \
      3A_S9/superexpressed.gtf \
      3C_S11/superexpressed.gtf \
      4A_S6/superexpressed.gtf \
      4C_S12/superexpressed.gtf \
      5A_S2/superexpressed.gtf \
      5C_S4/superexpressed.gtf \
      6A_S3/superexpressed.gtf \
      6C_S10/superexpressed.gtf
   if [ -e z1 ]; then
      LC_ALL=C sort -k 1,1 -k 4n,4 -k 3r,3 -k 5n,5 z1 > common_superexpressed.gtf
      rm z1
   fi
   LC_ALL=C sort -nk 1,1 z2 > coincidence_superexpressed.txt
   rm z2

   # There are 113 transcripts (321 exons) that are highly expressed in all samples.
   # Among them, 104 transcripts had already been identified as highly expressed in
   # all samples in the 2017-01-24 analysis.
fi

# Create a file called assemblies.txt that lists assembly files for each sample
if [ ! -e assemblies.txt ]; then
   find . -name transcripts.gtf > assemblies.txt
fi

# Create a simple merged transcriptome annotation.
if [ ! -e merged.gtf ]; then
   cuffmerge -p 24 \
             -g $ANNOTATION \
             -s $REFERENCE \
             -o ./ ./assemblies.txt 1> cuffmerge.log 2> cuffmerge.error
fi
