# Rice-PCA-SNPs



Welcome! This repo contains an analysis of SNP data from the RNA sequencing of many rice genomes. This workflow was adapted from the Genomics Laboratory course BIS 180L at UC Davis. This repo contains part of a wider analysis of the rice RNA seq data. More specifically, this repo contains:

* Making PCA and MDS plots from SNP data
* Assigning sub populations using fastStructure
* Exploring the relationship between fastStructure populations and PCA plots

The analyses conducted in this repo of the rice sequencing data will subsequently be used for a GWAS.


### File structure
---

This repo contains 3 directories: `input`, `output` and `scripts`. `input` and `output` contain input and output data to be fed into programs (like fastStructure), and `scripts` contains a .Rmd file with all code used in this analysis and the outputs of it, knitted to html format for viewing. In the future, I will also include .R scripts which will contain all code necessary to carry out the analyses.

### Workflow Overview
---

This analysis begins with loading data from a .csv file that contains processed SNP locations derived from a .vcf file, contaning sample information in the left column, with all subsequent columns as genome locations, and column entries as SNPs at those locations. Then:

1. Missing data was filled in with the average genotype
2. Principal components were computed, appended to the data
3. PCA and MDS plots were created, along with a plot showing the variance explained by each PC

The PCA plots suggested the existence of 3 or 4 sub populations in the 
