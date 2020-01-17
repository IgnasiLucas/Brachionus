# Differential expression in diapausing eggs 

This project contains some analysis performed by Dr. Eva Tarazona Castelblanque, for one chapter of her
PhD thesis, which is available here: http://roderic.uv.es/handle/10550/67254. It is a comparison of
gene expression levels in diapausing eggs from rotifer *Brachionus plicatilis*, from RNA-seq data.

In this repository, my aim is to improve the reproducibility of Eva's results and potentially contribute
some details to the analysis. Below, I will update the summaries of the results in reverse chronological
order. Each entry corresponds to one subfolder in the results folder, named after the day the analysis
started. The README.sh executable and documented files are meant to reproduce the results in each folder.

All the analyses are run in a Linux-64 platform. The file spec-file.txt is created by:

```bash
conda list --explicit > spec-file.txt
```

It can be used to create a conda environment with the same packages used originally to produce the
results. To re-create the conda environment that I call Brachionus do this:

```bash
conda create --name Brachionus --file spec-file.txt
```

Some folders have their own conda environment.

---------------------------------------------------------------------------------------------------

# 2020-01-14
Here I determine what gene ontology terms are associated with differential expression between
selective regimes, and between hatching conditions. Hatching condition did not affect many genes,
and a gene set enrichment analysis seems unnecessary in this case. However, for the sake of
completion, I report results for the following:

  * Genes differentially expressed between selective [regimes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-14/genes/regime/enrichment.html).
  * Genes differentially expressed between hatching [conditions](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-14/genes/hatching/enrichment.html).
  * Isoforms (or transcripts) differentially expressed between selective [regimes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-14/isoforms/regime/enrichment.html).
  * Isoforms (or transcripts) differentially expressed between hatching [conditions](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-14/isoforms/regime/enrichment.html).

# 2020-01-08
Properly modelling the experiment as a split-plot design is actually possible with the R
package `variancePartition`. Here I use it to analyse the proportion of variance of gene and
isoform expression explained by selective regime and by hatching treatment after properly
accounting for the random effect of the original population. I also use the `dream()` function
to test for significant fixed effects. This analysis is preferred to the previous one, and
will be used to generate the lists of genes and isoforms for functional analysis.

The reports for [genes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-08/genes/mixed.html)
and for [isoforms](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2020-01-08/isoforms/mixed.html)
are available.

# 2019-12-23
I examine in more detail the effect of the random environment in gene expression. Here I take
into account that not all populations respond in the same way to the selective regime. I fit a
new model which makes contrasts easier to interprete. Then I determine what and how many genes
are regulated in the random regime, relative to the regular one, individually in every population
under the random regime, and also in all three of them simultaneously. While several genes may
be regulated by the random environment in only one or two of the populations, most of the response
to this regime is common in all three populations. See the detailed report for either
[genes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-12-23/genes/randomness.html)
or [isoforms](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-12-23/isoforms/randomness.html).

# 2019-07-26
Use the R package topGO to run a functional enrichment analysis among the genes (or transcripts) that
are differentially expressed between selective regimes or hatching conditions. Find rendered versions
of the final reports here:

* [Genes differentially expressed between selective regimes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/genes/regime.html)
* [Genes differentially expressed between hatching conditions](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/genes/hatching.html)
* [Transcripts differentially expressed between selective regimes](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/transcripts/regime.html)
* [Transcripts differentially expressed between hatching conditions](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/transcripts/hatching.html)

# 2019-07-19
I use InterProScan to assign functional annotations to the proteins identified by TransDecoder among
the transcripts. Only 25% of the original 77728 transcripts (18% of the original genes) get any GO annotation.

# 2019-07-10
To start the functional annotation of transcripts, I use the gff2fasta tool of the cgat package to
extract the sequences of the 77728 transcripts. Then I use TransDecoder to identify the most promising
protein encoded within each transcript. TransDecoder predicts a protein in no more than 49663
transcripts.

# 2019-04-03
I use the package edgeR to identify genes differentially expressed between selective regimes and
hatching conditions. Find the rendered versions of the RMarkdown reports here:

* Using genes as units of expression:
  * [Preliminar analysis](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/genes/Preliminar.html)
  * [Effect of hatching condition](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/genes/Hatching.html)
  * [Interaction between hatching condition and selective regime](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/genes/Interactions.html)
  * [Effect of the selective regime](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/genes/Regime.html)
* Using isoforms (transcripts):
  * [Preliminar analysis](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/isoforms/Preliminar.html)
  * [Effect of hatching condition](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/isoforms/Hatching.html)
  * [Interaction between hatching condition and selective regime](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/isoforms/Interactions.html)
  * [Effect of the selective regime](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-04-03/isoforms/Regime.html)

# 2019-03-29
I quantify the expression using an alternative method. Here I use RSEM to obtain raw counts
of reads mapped to transcripts. This is the kind of input used by edgeR or DSeq2.

# 2019-03-21
Here I run cufflinks and cuffmerge to identify expressed transcripts. This is equivalent to what
Eva did on 2017-01-24. Here I made sure to include information on the type of library, since it is
stranded, with all reads expected on the reverse strand of the transcripts.

# 2019-03-19
Here I repeat the mapping of reads to the reference genome using Tophat. The main difference with
respect with the original mapping (Eva, 2016-12-19) is that I specify here the type of library. The
results are very much alike.

# 2019-02-01
Testing an alternative method of transcript quantification, for the sake of robustness of results. I run
kallisto, which is fast because it does not really map the reads to the transcripts, but performs a
pseudoalignment, based on k-mers. Kallisto did not work well: only ~30% of reads were pseudo-aligned.
The problem seems to be the fact that only about 67% of sequenced genomic fragments overlap with
predicted transcripts. Plus, the reference genome contains unknown sequences that could compromise the
pseudoalignment.

# 2017-12-29
Analysis of Gene Ontology terms.

# 2017-12-22
Lists of genes significantly differentially expressed (q < 0.05) with a minimum log(fold change) of 2.

# 2017-12-07
Incomplete analysis on the functional annotation of genes of interest.

# 2017-11-30
Comparisons of the lists of differentially expressed genes.

# 2017-10-27
Lists of differentially expressed genes.

# 2017-01-31
Estimation of expression levels and differential expression among conditions and regimes with Cuffdiff. 

# 2017-01-24
Assembly of transcripts with cufflinks and cuffmerge.

# 2016-12-19
Mapping of raw reads to reference genome.

| Sample | Reads    | Mapped   |  Rate   |
| ------ | --------:| --------:|:------ :|
| 1A_S8  | 58364853 | 46103494 | 78.99 % |
| 1C_S1  | 29222289 | 25433737 | 87.04 % |
| 2A_S7  | 55711941 | 49369286 | 88.62 % |
| 2C_S5  | 38563452 | 34283803 | 88.90 % |
| 3A_S9  | 55778100 | 48926714 | 87.72 % |
| 3C_S11 | 44386767 | 38807814 | 87.43 % |
| 4A_S6  | 59268714 | 51793879 | 87.39 % |
| 4C_S12 | 37364770 | 32119977 | 85.96 % |
| 5A_S2  | 39381051 | 30637627 | 77.80 % |
| 5C_S4  | 36745427 | 31601100 | 86.00 % |
| 6A_S3  | 55388301 | 48681614 | 87.89 % |
| 6C_S10 | 42730345 | 37107511 | 86.84 % |


# 2016-12-14
Quality control of RNA-seq reads. 

| Sample | MinLength | Average | MaxLength | NumSeqs  | Q20   | Q30   | MinQ | Average | MaxQ |
| ------ |:---------:|:-------:|:---------:|:--------:|:-----:|:-----:|:----:|:-------:|:----:|
| 1A_S8  |        20 |   74.39 |        75 | 58364853 | 96.10 | 94.31 |    2 |  34.44  |  36  |
| 1C_S1  |        21 |   74.39 |        75 | 29222289 | 96.12 | 94.39 |    2 |  34.46  |  36  |
| 2A_S7  |        20 |   74.43 |        75 | 55711941 | 96.12 | 94.36 |    2 |  34.45  |  36  |
| 2C_S5  |        22 |   74.50 |        75 | 38563452 | 96.19 | 94.46 |    2 |  34.47  |  36  |
| 3A_S9  |        20 |   74.46 |        75 | 55778100 | 96.12 | 94.37 |    2 |  34.45  |  36  |
| 3C_S11 |        21 |   74.48 |        75 | 44386767 | 96.11 | 94.37 |    2 |  34.45  |  36  |
| 4A_S6  |        21 |   74.45 |        75 | 59268714 | 96.11 | 94.34 |    2 |  34.45  |  36  |
| 4C_S12 |        20 |   74.48 |        75 | 37364770 | 95.99 | 94.17 |    2 |  34.41  |  36  |
| 5A_S2  |        23 |   74.47 |        75 | 39381051 | 95.64 | 93.79 |    2 |  34.33  |  36  |
| 5C_S4  |        20 |   74.48 |        75 | 36745427 | 95.92 | 94.09 |    2 |  34.39  |  36  |
| 6A_S3  |        20 |   74.48 |        75 | 55388301 | 96.14 | 94.40 |    2 |  34.46  |  36  |
| 6C_S10 |        22 |   74.49 |        75 | 42730345 | 95.99 | 94.23 |    2 |  34.42  |  36  |

