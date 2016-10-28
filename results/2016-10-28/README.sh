#!/bin/bash
#
# This is just a template, or example, with some initial tips.
#
# The README.sh file is meant to be executable. The very first line must
# be that one: "#!/bin/bash", to be interpreted by the bash shell. Once
# you run 'chmod u+x README.sh', you can call this file directly, to produce
# all the results, usually like this:
#
# ./README.sh 1> readme.log 2> readme.err &
#
# Directing the outputs is a good idea. In addition, the ampersand makes
# the script run in the background, so that you can keep using the terminal
# for other purposes, or even shut down the connexions while the server
# keeps working.
#
# It is a good idea to make the execution conditional, so that if the intended
# result is already present, there is no need to repeat that particular step
# of the pipeline. One way to make execution conditional is by checking the
# existence of the output file:

if [ ! -e z1.txt ]; then
   # Run a very long process here.
   echo something > z1.txt
fi

rm z1.txt

if [ ! -e summary.txt ]; then
   # Run the summary script to produce a summary table called 'summary.txt'.
   sleep 1
fi

if [ ! -e plot.png ]; then
   # I usually plot things with R, using a separate script, 'plot.R'.
   if [ -e plot.R ]; then
      R --no-save < plot.R
   else
      echo 'I cannot plot anything, because plot.R does not exist.'
   fi
fi

# Another useful thing is a list:

SAMPLE=(A1 A2 A3 A4 B1 B2 B3 B4)

#for i in 0 1 2 3 4 5 6 7; do
for i in `seq 0 7`; do
   echo "The name of the sample is ${SAMPLE[$i]}"
done

# Advanced stuff  :-D
# --------------
#
# The server has 64 processors, and plenty of RAM. We share it, and it is not
# a good idea to use all processors at the same time. But you can take advantage
# of this with simple parallelization tricks.

for i in 0 1 2 3 4 5 6 7; do
   # Run a long process with each sample, like this:
   if [ ! -e ${SAMPLE[$i]}.sam ] && [ -e reference.1.bt2 ]; then
      bowtie2 -x reference -U ${SAMPLE[$i]}.fastq -S ${SAMPLE[$i]}.sam 2> ${SAMPLE[$i]}.log &
   fi
done
wait

# If you parallelize like this, please add the word 'wait' after the loop, so
# that the script does not keep running before the previous processes finished.
#
# Also keep in mind that each process will require its own output, to prevent
# them from overwriting each other. This includes standard errors and standard
# outputs.
#
# It is a good idea to make each individual process conditional, so that you can
# delete and re-generate the results for a single sample without having to re-do
# the whole thing.
