# Differential expression in diapausing eggs 

This project contains some analysis performed by Dr. Eva Tarazona Castelblanque, for one chapter of her
PhD thesis, which is available here: http://roderic.uv.es/handle/10550/67254. It is a comparison of
gene expression levels in diapausing eggs from rotifer Brachionus plicatilis, from RNA-seq data.

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

---------------------------------------------------------------------------------------------------
# 2019-02-01
Testing an alternative method of transcript quantification, for the sake of robustness of results. I run
kallisto, which is fast because it does not really map the reads to the transcripts, but performs a
pseudoalignment, based on k-mers.

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

