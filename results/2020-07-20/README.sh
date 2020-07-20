#!/bin/bash
#
#				2020-07-20
#				==========
#
# Before submitting the manuscript, we need to prepare the data to make it available.
# I plan to archive the trimmed sequenced reads, rather than the very raw reads, because
# the trimmed reads are the starting point of my own analyses, and I consider them
# reliable enough.

DATADIR=../../data/trimmed

if [ ! -d ../../data/share ]; then mkdir ../../data/share; fi
if [ ! -e ../../data/share/trimmed.tar ]; then
   tar -cvf ../../data/share/trimmed.tar -C $DATADIR 1A_S8.trim.fastq.gz 1C_S1.trim.fastq.gz 2A_S7.trim.fastq.gz \
        2C_S5.trim.fastq.gz 3A_S9.trim.fastq.gz 3C_S11.trim.fastq.gz 4A_S6.trim.fastq.gz 4C_S12.trim.fastq.gz \
        5A_S2.trim.fastq.gz 5C_S4.trim.fastq.gz 6A_S3.trim.fastq.gz 6C_S10.trim.fastq.gz
   md5sum $DATADIR/* | sed 's,../../data/trimmed/,,' > trimmed_md5sum.txt
   tar -rf ../../data/share/trimmed.tar trimmed_md5sum.txt
fi

if [ ! -e ../../data/share/metadata.txt ]; then
   echo -e "#Sample\tRegime\tCondition\tPopulation"        > ../../data/share/metadata.txt
   echo -e "1A_S8\tunpredictable\tnon-forced_diapause\t1" >> ../../data/share/metadata.txt
   echo -e "1C_S1\tunpredictable\tforced_diapause\t1"     >> ../../data/share/metadata.txt
   echo -e "2A_S7\tunpredictable\tnon-forced_diapause\t2" >> ../../data/share/metadata.txt
   echo -e "2C_S5\tunpredictable\tforced_diapause\t2"     >> ../../data/share/metadata.txt
   echo -e "3A_S9\tpredictable\tnon-forced_diapause\t3"   >> ../../data/share/metadata.txt
   echo -e "3C_S11\tpredictable\tforced_diapause\t3"      >> ../../data/share/metadata.txt
   echo -e "4A_S6\tunpredictable\tnon-forced_diapause\t4" >> ../../data/share/metadata.txt
   echo -e "4C_S12\tunpredictable\tforced_diapause\t4"    >> ../../data/share/metadata.txt
   echo -e "5A_S2\tpredictable\tnon-forced_diapause\t5"   >> ../../data/share/metadata.txt
   echo -e "5C_S4\tpredictable\tforced_diapause\t5"       >> ../../data/share/metadata.txt
   echo -e "6A_S3\tpredictable\tnon-forced_diapause\t6"   >> ../../data/share/metadata.txt
   echo -e "6C_S10\tpredictable\tforced_diapause\t6"      >> ../../data/share/metadata.txt
fi

