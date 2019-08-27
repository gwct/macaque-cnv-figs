# Data, code, and figures for "Origins and long-term patterns of copy-number variation in rhesus macaques"

## Author
#### Gregg Thomas

## About

### This repository hosts the CNV calls and other data along with R scripts for analysis and reproducible figure generation.

## Citation

#### Forthcoming

## Contents

Each `fig` folder contains code to generate the figures in the paper as well as pre-rendered figures. Most figures can be generated in either color or greyscale if you change the `color_plots` variable near the top of each script.

Here is a description of each folder:

#### 1. `counts`: This folder contains a script that just reads and reports CNV counts in a log file. 
#### 2. `data`: All the data used to generate figures in the paper

* `brandler-cafe-genes.csv`: Ensembl IDs for human genes that overlap with CNV calls from Brandler et al.'s data.
* `brandler-cnvs.zip`: CNV calls from Brandler et al. (2016) in long format with lots of info. This is the main file used for analysis of human CNVs. This file has been zipped because the original was over github's 100MB file size limit. These calls are unfiltered from their supplemental data S1 (their filtering is replicated in the scripts).
* `brandler-denovo.csv`: The *de novo* CNVs called by Brandler et al. (2016) for humans.
* `brandler-samples.csv`: Individual and family identifiers for the Brandler et al, (2016) human data.
* `cafe-traits.csv`: Statistics for the gene family analysis of 17 mammals.
* `cafe-tree.tre`: The newick formatted phylogeny for the 17 mammals used in the CAFE analysis.
* `macaque-cafe-genes-filtered.csv`: Ensembl IDs for macaque genes that overlap with CNV calls after filtering.
* `macaque-cnv-chrome-counts.csv`: The number and type of CNVs called for each macaque chromosome.
* macaque-cnvs.csv: The macaque CNV calls in long format. Unfiltered (filtering is done in the scripts).
* `macaque-peds.pdf`: Pedigree structures and individual IDs for the 32 rhesus macaques sampled.

#### 3. `fig1`: Proportions of CNV types and bases for both humans and macaques as well as a general pedigree structure for our macaque sample.
#### 4. `fig2`: Locations of macaque CNVs along chromosomes.
#### 5. `fig3`: Length distributions of human and macaque CNVs.
#### 6. `fig4`: *de novo* CNVs in humans and macaques.
#### 7. `fig5`: Phylogeny of 17 mammals used in gene family analysis and proportion of gene count changes that are fixed and segregating in macaques.
#### 8. `figS1`: Correlation between chromosome length and number of CNVs.
#### 9. `lib`: Supplementary scripts that read and filter CNV data.
#### 10. `other-figs`: A couple of other figures I generated but didn't include in the paper.
* `figX1`: No correlation between CNV length and number of genes overlapped.
* `figX2`: CNV changes vs. gene gains and losses

## Reference

The human CNV data was obtained from:

Brandler WM, Antaki D, Gujral M, Noor A, Rosanio G, Chapman TR, Barrera DJ, Lin GN, Malhorta D, Watts AC, Wong LC, Estabillo JA, Gadomski TE, Hong O, Fajardo KVF, Bhandari A, Owen R, Baughn M, Yuan J, Solomon T, Moyzis AG, Maile MS, Sanders SJ, Reiner GE, Vaux KK, Strom CM, Zhang K, Muotri AR, Akshoomoff N, Leal SM, Pierce K, Courchesne E, Iakoucheva LM, Corsello C, Sebat J. 2016. Frequency and complexity of de novo structural mutation in autism. *Am J Hum Genet*. 98(4):667-79. doi: https://doi.org/10.1016/j.ajhg.2016.02.018