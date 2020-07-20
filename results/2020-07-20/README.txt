
Gene expression in rotifer diapausing eggs in response to divergent environmental predictability regimes

       2020. Eva Tarazona, J. Ignacio Lucas-Lledó, María José Carmona & Eduardo M. García-Roger

eduardo.garcia@uv.es

----------------------------------------------------------------------------------------------

The present README.txt file accompanies two files: trimmed.tar and metadata.txt.

The trimmed.tar file is an archive containing 12 compressed fastq files,
and the trimmed_md5sum.txt file. To extract them, run:

   tar -xf trimmed.tar

The trimmed_md5sum.txt can be used to check the integrity of the fastq.gz
files, like this:

   md5sum -c trimmed_md5sum.txt

The fastq.gz files are named by the sample name. Samples are described
in the metadata.txt file. The fastq.gz files contain single end reads from
RNA-seq experiments, after trimming and quality filtering.

The open notebook at https://github.com/ignasilucas/Brachionus describes the
analysis performed with this data by the authors. Feel free to use the data under
the terms of the license. You can address questions about the data set to
eduardo.garcia@uv.es.


