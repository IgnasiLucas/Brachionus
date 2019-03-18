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
KMER_SIZE=23

if [ ! -e summary.txt ]; then
   if [ ! -e chromosomes.txt ]; then
      gawk '($3 == "contig"){print $1 "\t" $5}' $ANNOTATION | LC_ALL=C sort -nrk 2,2 > chromosomes.txt
   fi
   if [ ! -e index ]; then
      if [ ! -e transcriptome.fa ]; then
         if [ ! -e mRNA.gff ]; then
            # I remove duplicated mRNA, and order them.
            grep -P "\tmRNA\t" $ANNOTATION | \
            gawk '($1 ":" $4 "-" $5 != PREV){print}{PREV = $1 ":" $4 "-" $5}' | \
            bedtools sort -g chromosomes.txt -i - > mRNA.gff
         fi
         bedtools getfasta -fi $REFERENCE -bed mRNA.gff -fo transcriptome.fa
      fi
      kallisto index -i index --make-unique --kmer-size $KMER_SIZE transcriptome.fa
   fi
   for i in $(seq 0 11); do
      if [ ! -d ${SAMPLES[i]} ]; then mkdir ${SAMPLES[i]}; fi
      if [ ! -e ${SAMPLES[i]}/abundance.h5 ]; then
         kallisto quant -i index -o ${SAMPLES[i]} --single -l 200 -s 50 \
            --bias --threads 1 --single-overhang --rf-stranded --genomebam \
            --gtf mRNA.gff --chromosomes chromosomes.txt $FASTQDIR/${SAMPLES[$i]}.trim.fastq.gz
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

# Only about 30% of reads get pseudoaligned. This contrasts with the high (~80%) mapping success
# obtained originally by Eva (2016-12-19). Some reads may fail to map here because I enforce the stranded
# (pseudo)mapping. I wonder if using the largest k-mers (31-mers) may also reduce the pseudomapping
# success. I will run it again with shorter k-mers.
#
# After trying a couple different values of k-mer size, I realize the main reason for the low
# (pseudo)mapping success here is that kallisto is mapping reads only to the known transcripts,
# not anywhere else in the genome, which may not be the most sensitive thing to do in this
# case, given that structural annotation of Brachionus plicatilis genome is based on predictions
# by the program 'maker'. Many real transcripts may not be found among the predicted ones. It is actually
# a good question and easier to answer: how many (and how much) of the predicted transcripts overlap
# with sequenced fragments? And the other way around: how many of the apparently expressed genome
# overlap predicted transcripts? During alignment of the reads, Tophat was already aware of the
# reference set of predicted transcripts, and used it as a preferred target for mapping, before
# mapping to the rest of the genome. Unfortunately, Tophat did not report how many reads mapped to
# the predicted genes. Then, cufflinks was also informed of the predicted annotation,
# and in the output of cufflinks all reference transcripts are supposed to be included.

if [ ! -e expressed_predicted.bed ]; then
   if [ ! -e expressed.bed ]; then
      if [ ! -e expressed.bam ]; then
         samtools merge -r expressed.bam ../2016-12-19/1A_S8/accepted_hits.bam \
                                         ../2016-12-19/1C_S1/accepted_hits.bam \
                                         ../2016-12-19/2A_S7/accepted_hits.bam \
                                         ../2016-12-19/2C_S5/accepted_hits.bam \
                                         ../2016-12-19/3A_S9/accepted_hits.bam \
                                         ../2016-12-19/3C_S11/accepted_hits.bam \
                                         ../2016-12-19/4A_S6/accepted_hits.bam \
                                         ../2016-12-19/4C_S12/accepted_hits.bam \
                                         ../2016-12-19/5A_S2/accepted_hits.bam \
                                         ../2016-12-19/5C_S4/accepted_hits.bam \
                                         ../2016-12-19/6A_S3/accepted_hits.bam \
                                         ../2016-12-19/6C_S10/accepted_hits.bam
      fi
      bedtools merge -i expressed.bam -s | bedtools sort -g chromosomes.txt -i - > expressed.bed
   fi
   # The mRNA.gff are the predicted exons.
   bedtools intersect -a expressed.bed -b mRNA.gff -c -sorted -g chromosomes.txt > expressed_predicted.bed
fi

if [ ! -e predicted_expressed.gff ]; then
   bedtools intersect -a mRNA.gff -b expressed.bed -c -sorted -g chromosomes.txt > predicted_expressed.gff
fi

if [ ! -e summary2.txt ]; then
   NUM_EXPRESSED_PREDICTED=$( gawk '($4 > 0){N++}END{print N}' expressed_predicted.bed )
   NUM_EXPRESSED=$( cat expressed.bed | wc -l )
   NUM_PREDICTED_EXPRESSED=$( gawk '($NF > 0){N++}END{print N}' predicted_expressed.gff )
   NUM_PREDICTED=$( cat mRNA.gff | wc -l )
   echo -e "Number of expressed exons:               $NUM_EXPRESSED"            > summary2.txt
   echo -e "Number of expressed and predicted exons: $NUM_EXPRESSED_PREDICTED" >> summary2.txt
   echo -e "Number of predicted genes:               $NUM_PREDICTED"           >> summary2.txt
   echo -e "Number of predicted and expressed genes: $NUM_PREDICTED_EXPRESSED" >> summary2.txt
fi

# About 67% of expressed exons had been predicted by maker. And 95% of predicted genes have at least
# one read mapping to them:
#
#     Number of expressed exons:               277504
#     Number of expressed and predicted exons: 188096
#     Number of predicted genes:                54870
#     Number of predicted and expressed genes:  52132
#
# There are 89408 expressed exons that are not included in the annotation.

if [ ! -e not_annotated_expressed.bed ]; then
   gawk '($4 == 0){print}' expressed_predicted.bed > not_annotated_expressed.bed
fi
