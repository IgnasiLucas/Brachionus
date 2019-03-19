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
MERGED=../2017-01-24/cuffmerge/merged.gtf
FASTQDIR=../../data/trimmed
SAMPLES=( 1A_S8 1C_S1 2A_S7 2C_S5 3A_S9 3C_S11 4A_S6 4C_S12 5A_S2 5C_S4 6A_S3 6C_S10 )
KMER_SIZE=23

if [ ! -e summary.txt ]; then
   if [ ! -e chromosomes.txt ]; then
      gawk '($3 == "contig"){print $1 "\t" $5}' $ANNOTATION | LC_ALL=C sort -nrk 2,2 > chromosomes.txt
   fi
   if [ ! -e index ]; then
      if [ ! -e transcriptome.fa ]; then
         if [ ! -e exons.gtf ]; then
            bedtools sort -g chromosomes.txt -i $MERGED > exons.gtf
         fi
         bedtools getfasta -fi $REFERENCE -bed exons.gtf -fo transcriptome.fa
      fi
      kallisto index -i index --make-unique --kmer-size $KMER_SIZE transcriptome.fa
   fi
   for i in $(seq 0 11); do
      if [ ! -d ${SAMPLES[i]} ]; then mkdir ${SAMPLES[i]}; fi
      if [ ! -e ${SAMPLES[i]}/abundance.h5 ]; then
         kallisto quant -i index -o ${SAMPLES[i]} --single -l 200 -s 80 \
            --bias --threads 1 --single-overhang --rf-stranded --genomebam \
            --gtf exons.gtf --chromosomes chromosomes.txt $FASTQDIR/${SAMPLES[$i]}.trim.fastq.gz
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

# CONCLUSION
# ==========
#
# I definitely cannot use kallisto, much to my desappointment. First, I tried using the predicted
# transcripts. That produced very low mapping success (around 30% of reads). One main reason for
# this must be the fact that many reads do not map to predicted transcripts. Only around 67% of
# expressed sequenced fragments overlap predicted transcripts. There may be other reasons, such
# as the incompleteness of the reference genome. It seems that Kallisto substitutes Ns in the
# reference for pseudo-random nucleotides. I don't know what effect that has on the pseudoallignment,
# but I imagine it could precent the pseudo-mapping of some reads.
#
# To alleviate the first problem, I used as reference the merged transcriptome assembly obtained
# with Cufflinks and Cuffmerge by Eva in 2017-01-24. However, the mapping success hardly improved.
# One obvious problem with this approach is that the reference transcriptome built like that is
# composed of exons, not full transcripts. Some exons are very short. Reads encompassing more than
# one exon may not get properly mapped, and the whole quantification must be screwed anyways.
#
# I erase the results to save space.
