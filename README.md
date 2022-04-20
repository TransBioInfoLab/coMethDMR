## coMethDMR: Accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies 
Gomez L, Odom GJ, Young JI, Martin ER, Liu L, Chen X, Griswold AJ, Gao Z, Zhang L, Wang L (2019) Nucleic Acids Research, gkz590, https://doi.org/10.1093/nar/gkz590



## Description
coMethDMR is an R package that identifies genomic regions associated with continuous phenotypes by optimally leverages covariations 
among CpGs within predefined genomic regions. Instead of testing all CpGs within a genomic region, coMethDMR carries out an additional 
step that selects comethylated sub-regions first without using any outcome information. Next, coMethDMR tests association between 
methylation within the sub-region and continuous phenotype using a random coefficient mixed effects model, which models both variations 
between CpG sites within the region and differential methylation simultaneously.



## Quick Overview / Analysis Steps
Assuming you have completed all pre-processing and normalization procedures, here are the steps to analyse your Illumina EPIC or 450k DNA methylation data.

1. **Your Data**:
    + Ensure your methylation data is loaded into R's Global Environment as a numeric matrix in *probe by sample* form: with probe IDs as your row names and sample IDs as the column names
    + Your phenotype / clinical data should be a data frame with a column called `Sample` (spelled exactly with an upper-case "S"); these sample IDs should match the column names of the methylation data
2. **Pre-Defined Regions**: Load the list of pre-calculated regions of "contiguous" CpGs which matches your Illumina data type. We have pre-calculated some of these lists of regions. We used the `CloseBySingleRegion()` function with `maxGap = 200` (genomic locations within 200 base pairs are placed in the same cluster) and `minCpGs = 3` (we need at least 3 CpGs to retain the location). These data files are:
    + Genic regions, 450k array: `extdata/450k_Gene_3_200.rds`. Load this via `system.file("extdata", "450k_Gene_3_200.rds", package = "coMethDMR", mustWork = TRUE)`
    + Inter-genic regions, 450k array: `extdata/450k_InterGene_3_200.rds`.Load this via `system.file("extdata", "450k_InterGene_3_200.rds", package = "coMethDMR", mustWork = TRUE)`
    + Genic regions, EPIC array: download the supplemental data file from <https://github.com/TransBioInfoLab/coMethDMR_data/blob/main/data/EPIC_10b4_Gene_3_200.rds>
    + Inter-genic regions, EPIC array: download the supplemental data file from <https://github.com/TransBioInfoLab/coMethDMR_data/blob/main/data/EPIC_10b4_InterGene_3_200.rds>
3. **Adjust Methylation for Covariates** with the `GetResiduals()` function; your methylation values may be confounded by clinical variables unrelated to your treatment, such as sex, age, or even [the square of age](https://www.nature.com/articles/s41598-021-88504-0)
4. **Regions of Co-Methylation**



## Installation


### Bioconductor Version
The `coMethDMR::` package has been accepted to the [Bioconductor](https://bioconductor.org/) repository of R packages. It will be included in version 3.15 (April 2022 release). To install [this version](https://www.bioconductor.org/packages/devel/bioc/html/coMethDMR.html), please use the following code:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("coMethDMR")
```


### Development Version
The development version of `coMethDMR::` can also be installed from this GitHub repository by

```{r eval=FALSE, message=FALSE, warning=FALSE, results='hide'}
library(devtools)
install_github("TransBioInfoLab/coMethDMR")
```

Please note that using compiled code from GitHub may require your computer to have additional software ([Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) for Windows or [Xcode](https://developer.apple.com/xcode/) for Mac). Also note that installing this development version may result in some errors. We have outlined potential troubleshooting steps below.

#### Install Errors: Cache
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

If so, please fix this by running `ExperimentHub::ExperimentHub()` first (and type `yes` if you receive a prompt to create a local cache for your data), then re-installing the package. Please see this white paper for more information: <https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html>.

#### Install Errors: `.onLoad()` Failure
You may also get this error during installation:
```
Error: package or namespace load failed for 'coMethDMR':
 .onLoad failed in loadNamespace() for 'coMethDMR', details:
  call: NULL
  error: $ operator is invalid for atomic vectors
```

This error is caused by a version mismatch issue for the `sesameData::` (<https://bioconductor.org/packages/sesameData/>) package. We require `sesameData::` version 1.12 or higher. To fix this, you will need Biocdonductor version 3.14 or later. The following code will assist here:
```
BiocManager::install(version = "3.14")
BiocManager::install("sesameData")
```

After successfully executing the above installation, you should be able to install `coMethDMR::` from GitHub like normal.


### Loading the Package
After installation, the coMethDMR package can be loaded into R using:

```{r eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(coMethDMR)
```



## Manual

The reference manual for coMethDMR can be downloaded from old repository: <https://github.com/TransBioInfoLab/coMethDMR_old/tree/master/docs/>. The reference manual is [coMethDMR_0.0.0.9001.pdf](https://github.com/TransBioInfoLab/coMethDMR_old/blob/master/docs/coMethDMR_0.0.0.9001.pdf). Two vignettes are available in the same directory: [1_Introduction_coMethDMR_10-9-2019.pdf](https://github.com/TransBioInfoLab/coMethDMR_old/blob/master/docs/1_Introduction_coMethDMR_10-9-2019.pdf) and [2_BiocParallel_for_coMethDMR_geneBasedPipeline.pdf](https://github.com/TransBioInfoLab/coMethDMR_old/blob/master/docs/2_BiocParallel_for_coMethDMR_geneBasedPipeline.pdf).



## Frequently Asked Questions

1. There are two main steps in coMethDMR: (1) identifying comethylatyed clusters (2) testing methylation levels in those comethylated clusters against a phenotype".
In step 1 and 2, should we use beta value or M values for CoMethDMR?

Answer: In step (1), using M values and beta values produce similar results. See Supplementary Table 2 Comparison of using beta values or M-values for identifying co-methylated regions in first step of coMethDMR pipeline at optimal rdrop parameter value of the coMethDMR paper. 

In step (2), M-values should be used because it has better statistical properties. See Du et al. (2010) Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. 



## Development History
Our development history is at https://github.com/TransBioInfoLab/coMethDMR_old

