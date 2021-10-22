## coMethDMR: Accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies 
Gomez L, Odom GJ, Young JI, Martin ER, Liu L, Chen X, Griswold AJ, Gao Z, Zhang L, Wang L (2019) Nucleic Acids Research, gkz590, https://doi.org/10.1093/nar/gkz590

## Description
coMethDMR is an R package that identifies genomic regions associated with continuous phenotypes by optimally leverages covariations 
among CpGs within predefined genomic regions. Instead of testing all CpGs within a genomic region, coMethDMR carries out an additional 
step that selects comethylated sub-regions first without using any outcome information. Next, coMethDMR tests association between 
methylation within the sub-region and continuous phenotype using a random coefficient mixed effects model, which models both variations 
between CpG sites within the region and differential methylation simultaneously.

## Installation

The latest version can be installed by

```{r eval=FALSE, message=FALSE, warning=FALSE, results='hide'}
library(devtools)
install_github("TransBioInfoLab/coMethDMR")
```
After installation, the coMethDMR package can be loaded into R using:

```{r eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(coMethDMR)
```

### Install Errors
You may get the following error during installation:
```
Error: package or namespace load failed for 'coMethDMR':
 .onLoad failed in loadNamespace() for 'coMethDMR', details:
  call: .updateHubDB(hub_bfc, .class, url, proxy, localHub)
  error: Invalid Cache: sqlite file
  Hub has not been added to cache
  Run again with 'localHub=FALSE'
Error: loading failed
```

If so, please fix this by running `ExperimentHub::ExperimentHub()` first (and type `yes` if you receive a prompt to create a local cache for your data), then re-installing the package.


## Manual

The reference manual for coMethDMR can be downloaded from old repository https://github.com/TransBioInfoLab/coMethDMR_old/tree/master/docs/coMethDMR_0.0.0.9001.pdf. Two vignettes are available in the same directory: "1_Introduction_coMethDMR_10-9-2019.pdf" and "2_BiocParallel_for_coMethDMR_geneBasedPipeline.pdf"

## Frequently Asked Questions

1. There are two main steps in coMethDMR: (1) identifying comethylatyed clusters (2) testing methylation levels in those comethylated clusters against a phenotype".
In step 1 and 2, should we use beta value or M values for CoMethDMR?

Answer: In step (1), using M values and beta values produce similar results. See Supplementary Table 2 Comparison of using beta values or M-values for identifying co-methylated regions in first step of coMethDMR pipeline at optimal rdrop parameter value of the coMethDMR paper. 

In step (2), M-values should be used because it has better statistical properties. See Du et al. (2010) Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. 


## Development History
Our development history is at https://github.com/TransBioInfoLab/coMethDMR_old

