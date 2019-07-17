#!/bin/bash
#
#				2019-07-10
#				==========
#
# Once the genes and transcripts differentially expressed between selective
# regimes are identified (2019-04-03), I need to run the functional analysis.
# Genes in the B. plicatilis genome do not seem to have any functional annotation.
# I believe that Eva used Blast2GO, but I'd like to try something else. I note
# that the complete set of transcripts, including those discovered along this
# project, are in 2019-03-29/z1.gtf. This file does not have CDS information, but
# it includes all the genes and transcripts used in the differential gene expression
# analysis, properly identified with the names given by cuffmerge.
#
# The cleanest way to run the functional analysis is to start from a fasta file
# with all transcripts. While bedtools allows to extract sequences from gff files,
# it does not join exons from the same transcript together (unless using a BED 12
# file). I have looked for "gff2fasta.py" scripts, and I learned about cgat scripts,
# which can be installed with conda. They are not compatible with the current
# environment. Thus, to run this folder you need to activate the cgat environment,
# which is saved here as env-cgat.yaml, or spec-file.txt.

GFF="../2019-03-29/z1.gtf"
REF="../../data/reference.fa"

# gff2fasta outputs a header in the fasta file, which I want to re-direct to a log file.
# It can also include transcript attributes found in the gtf file in the name of the
# sequences in the fasta file. I want to process the output of gff2fasta in two different
# and simultaneous ways: redirect the header to a log, and remove unnecessary attributes
# from the names. I learned that I can do this with a named pipe and the tee command.

if [ ! -e transcripts.fa ]; then
   mkfifo pipe
   # Here, I make grep read the pipe in the background, before the pipe carries anything.
   cat pipe | grep "^#" > transcripts.log &
   # Now, I send the output of gff2fasta both to the pipe and to the other grep command.
   cgat gff2fasta --genome-file $REF \
                  --merge-adjacent \
                  --is-gtf \
                  --header-attributes < $GFF | \
   tee pipe | \
   # I expect only the name lines of the fasta file to have spaces, and I expect the
   # first and the forth fields to be transcript and gene ids, respectively.
   grep -v "^#" | cut -d " " -f 1,4 > transcripts.fa
   rm pipe
fi

if [ ! -e transdecoder/longest_orfs.pep ]; then
   TransDecoder.LongOrfs -t transcripts.fa -S -O transdecoder
fi

if [ ! -e blastp.outfmt6 ]; then
   if [ ! -e swissprot.00.phr ]; then
      update_blastdb.pl --decompress swissprot
   fi
   blastp -query transdecoder/longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
fi
