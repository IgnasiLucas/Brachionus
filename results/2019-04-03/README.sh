#!/bin/bash
#
#				2019-04-03
#				==========
#
# Here I will run edgeR to identify differentially expressed genes. The samples
# consist of 3 replicates of the 4 combinations of two two-level treatments:
# the selective regime (predictable or unpredictable) and the hatching time
# (either after a forced 28 days diapause, or right away). The distribution of
# samples among treatments was the following:
#
#  +------------------------------------+------------------------------------+
#  |    Predictable selective regime    |  Unpredictable selective regime    |
#  +------------------------------------+------------------------------------+
#  | Forced diapause | Without diapause | Forced diapause | Without diapause |
#  +-----------------+------------------+-----------------+------------------+
#  |     3C_S11      |      3A_S9       |     1C_S1       |      1A_S8       |
#  |     5C_S4       |      5A_S2       |     2C_S5       |      2A_S7       |
#  |     6C_S10      |      6A_S3       |     4C_S12      |      4A_S6       |
#  +-----------------+------------------+-----------------+------------------+
#
# To put it another way:
#
#  +--------+------------------+------------------+
#  | Sample | Selective regime | Forced diapause? |
#  +--------+------------------+------------------+
#  | 1A_S8  | Unpredictable    | No               |
#  | 1C_S1  | Unpredictable    | Yes              |
#  | 2A_S7  | Unpredictable    | No               |
#  | 2C_S5  | Unpredictable    | Yes              |
#  | 3A_S9  | Predictable      | No               |
#  | 3C_S11 | Predictable      | Yes              |
#  | 4A_S6  | Unpredictable    | No               |
#  | 4C_S12 | Unpredictable    | Yes              |
#  | 5A_S2  | Predictable      | No               |
#  | 5C_S4  | Predictable      | Yes              |
#  | 6A_S3  | Predictable      | No               |
#  | 6C_S10 | Predictable      | Yes              |
#  +--------+------------------+------------------+
#
# Note that samples are paired between diapause treatments. After the selection
# experiment, "diapausing eggs obtained for each laboratory population were divided
# into two groups to be subjected to two different diapause conditions: (1) non-forced
# diapause (NFD) and (2) forced diapause (FD)." Thus, when testing for the effect
# of diapause (or hatching) condition, it would not be efficient to compare the
# average level of expression among samples forced to hatch with the average level
# among samples forced to diapause.
#
# So, I guess, the rigth description of the data includes an additional column:
#
#  +--------+------------------+------------------+-------+
#  | Sample | Selective regime | Forced diapause? | Block |
#  +--------+------------------+------------------+-------+
#  | 1A_S8  | Unpredictable    | No               |   1   |
#  | 1C_S1  | Unpredictable    | Yes              |   1   |
#  | 2A_S7  | Unpredictable    | No               |   2   |
#  | 2C_S5  | Unpredictable    | Yes              |   2   |
#  | 3A_S9  | Predictable      | No               |   3   |
#  | 3C_S11 | Predictable      | Yes              |   3   |
#  | 4A_S6  | Unpredictable    | No               |   4   |
#  | 4C_S12 | Unpredictable    | Yes              |   4   |
#  | 5A_S2  | Predictable      | No               |   5   |
#  | 5C_S4  | Predictable      | Yes              |   5   |
#  | 6A_S3  | Predictable      | No               |   6   |
#  | 6C_S10 | Predictable      | Yes              |   6   |
#  +--------+------------------+------------------+-------+

