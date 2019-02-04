#!/bin/bash
#
#				2019-02-01
#				==========
#
# There are several methods available to analyse RNA-seq data. It is clear that not a single
# method is best for all data sets. Here, I use an innovative method by the Pachter lab, called
# kallisto to quantify transcript abundances (Bray et al. 2016, Nat. Biotech. 34:525-527). I
# instruct kallisto to map reads on the reverse strand of transcripts only (--rf-stranded),
# because the TrueSeq Illumnia technology produces stranded reads that should map on the reverse
# strand.

REFERENCE=../../data/reference.fa
ANNOTATION=../../data/annotation.gff
FASTQDIR=../../data/trimmed
SAMPLES=( 1A_S8 1C_S1 2A_S7 2C_S5 3A_S9 3C_S11 4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )

if [ ! -e summary.txt ]; then
   if [ ! -e index ]; then
      if [ ! -e transcriptome.fa ]; then
         if [ ! -e mRNA.gff ]; then
            grep -P "\tmRNA\t" $ANNOTATION > mRNA.gff
         fi
         bedtools getfasta -fi $REFERENCE -bed mRNA.gff -fo transcriptome.fa
      fi
      kallisto index -i index --make-unique transcriptome.fa
   fi
   if [ ! -e chromosomes.txt ]; then
      gawk '($3 == "contig"){print $1 "\t" $5}' $ANNOTATION > chromosomes.txt
   fi
   for i in $(seq 0 11); do
      if [ ! -d ${SAMPLES[i]} ]; then mkdir ${SAMPLES[i]}; fi
      if [ ! -e ${SAMPLES[i]}/abundance.h5 ]; then
         kallisto quant -i index -o ${SAMPLES[i]} --single -l 200 -s 50 \
            --bias --threads 1 --single-overhang --rf-stranded --genomebam \
            --gtf $ANNOTATION --chromosomes chromosomes.txt $FASTQDIR/${SAMPLES[$i]}.trim.fastq.gz
      fi
   done
   gawk 'BEGIN{
      print "sample\tNum_reads\tPseudoaligned\tPercent"
   }(/will process file/){
      FASTQ = substr($6,20,2)
   }(/reads pseudoaligned/){
      gsub(/,/,"")
      printf("%s\t% 8i\t% 9i\t%.2f\n", FASTQ, $3, $5, 100 * $5 / $3)
   }' readme.err > summary.txt
fi

# Only about 30% of reads get pseudoaligned. This contrasts with the large mapping success
# obtained originally by Eva. Some reads may fail to map here because I enforce the stranded
# (pseudo)mapping. I wonder if using the largest k-mers (31-mers) can reduce the pseudomapping
# success. I will run it again with shorter k-mers.
