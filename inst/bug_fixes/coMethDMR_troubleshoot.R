library(tidyverse)
betas <- read_delim("~/Desktop/betas.txt",delim=" ")
pheno <- read_delim("~/Desktop/pheno.txt",delim=" ")
cometh <- readRDS("~/Downloads/coMeth_16.rds")
cometh[["chr16:2961485-2962051"]]
filtered <- betas %>% 
  filter(CpG %in% cometh[["chr16:2961485-2962051"]]) %>% 
  as.data.frame()

#These are the other comethylated regions I included in the betas file.#
cometh[["chr16:85645063-85645310"]]
cometh[["chr16:18801545-18801821"]]


row.names(filtered) <- filtered$CpG
filtered$rrow <- NULL
filtered$CpG <- NULL
filtered[,1:5]
names(pheno)[2] <- "Sample"

#Vector of covariates
fixed_cov=c("age","BMI","sex","PC1","PC2","PC3","PC4","CD8T","CD4T","NK","Bcell","Mono","batch")

library(coMethDMR)
lmmTest(beta=filtered,pheno,contPheno_char="LVMindex",covariates_char=fixed_cov,modelType="randCoef",arrayType="EPIC")



#This was my initial code to get the comethylated regions. No issues with this function.
gene_ls <- readRDS(system.file("extdata","EPIC_Gene_3_200.rds",package = 'coMethDMR',mustWork = TRUE))
coMeth <- CoMethAllRegions(dnam = x, betaToM=T, CpGs_ls = gene_ls, arrayType = "EPIC", returnAllCpGs=F, output = "CpGs")
saveRDS(coMeth,out.file)
#And then lmmTestAllRegions where I didn't get any beta estimates, std err, etc.
fit <- lmmTestAllRegions(betas=betas,region_ls=cometh,pheno,contPheno_char="LVMindex",covariates_char=fixed_cov,modelType="randCoef",arrayType="EPIC")
out <- AnnotateResults(lmmRes_df=fit,arrayType="EPIC")
