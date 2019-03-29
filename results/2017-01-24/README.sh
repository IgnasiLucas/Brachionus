#!/bin/bash

#Assemble expressed genes and transcripts

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
                 -N \
		 -q \
		 -g $ANNOTATION/maker2.all.gff \
		 -I 15000 \
		 --min-intron-length 10 \
                 -o ./$dirname.cufflinks $TOPHAT/$dirname/accepted_hits.bam 1> $dirname.cufflinks/log 2> $dirname.cufflinks/error
   fi
   if [ ! -e $dirname.cufflinks/fpkm_hist.png ]; then
      gawk '{gsub(/;|\"/,""); for (i=9; i<=NF; i++) {if ($i == "FPKM") print $(i+1)}}' $dirname.cufflinks/transcripts.gtf > z1
      Rscript plot.R $dirname.cufflinks/fpkm_hist.png z1
      rm z1
   fi
   # After looking at the histograms of relative expression level, I notice that all samples have
   # a kind of long tail of a few highly expressed genes. I arbitrarily establish the threshold of
   # highly expressed genes as those with PFKM > 1000. Below, I extract a list of those transcripts
   # and exons for every sample, so that we can compare them among samples. If there was a common
   # set of highly expressed genes in all samples, I would consider removing them from analysis with
   # the --mask-file option to cufflinks.
   if [ ! -e $dirname.cufflinks/superexpressed.gtf ]; then
      gawk -v FS="\t" '{
         gsub(/;|\"/,"")
         split($9,A,/ /)
         for (i=1; i<=length(A); i++) {
            if ((A[i] == "FPKM") && (A[i+1] + 0 > 1000.0)) print $0
         }
      }' $dirname.cufflinks/transcripts.gtf > $dirname.cufflinks/superexpressed.gtf
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
   }' 1A_S8.cufflinks/superexpressed.gtf \
      1C_S1.cufflinks/superexpressed.gtf \
      2A_S7.cufflinks/superexpressed.gtf \
      2C_S5.cufflinks/superexpressed.gtf \
      3A_S9.cufflinks/superexpressed.gtf \
      3C_S11.cufflinks/superexpressed.gtf \
      4A_S6.cufflinks/superexpressed.gtf \
      4C_S12.cufflinks/superexpressed.gtf \
      5A_S2.cufflinks/superexpressed.gtf \
      5C_S4.cufflinks/superexpressed.gtf \
      6A_S3.cufflinks/superexpressed.gtf \
      6C_S10.cufflinks/superexpressed.gtf
   if [ -e z1 ]; then
      LC_ALL=C sort -k 1,1 -k 4n,4 -k 3r,3 -k 5n,5 z1 > common_superexpressed.gtf
      rm z1
   fi
   LC_ALL=C sort -nk 1,1 z2 > coincidence_superexpressed.txt
   rm z2

   # There are 112 transcripts (330 exons) that are highly expressed in all samples.
fi

# Create a file called assemblies.txt that list assembly file for each sample

	if [ ! -e assemblies.txt ]; then
           find . -path "*/transcripts.gtf" > assemblies.txt #Para que liste todos los archivos
        fi

# Cuffmerge-Create a simple merged transcriptome annotation #input seria */transcripts.gtf
   if [ ! -d cuffmerge ]; then mkdir cuffmerge; fi
   if [ ! -e cuffmerge/merged.gtf ]; then
	cuffmerge -p 24 \
		  -g $ANNOTATION/maker2.all.gff \
                  -s $REFDIR/Brachionus_plicatilis_scaffold_min500.fasta \
		  -o cuffmerge ./assemblies.txt 1> cuffmerge.log 2> cuffmerge.error
   fi

#Cudiff-quantifies expresSion

