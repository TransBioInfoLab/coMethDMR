# Test Issue
# Gabriel Odom and Fernanda Veitzman
# 2021-04-27

# We received an issue on GitHub:
#   https://github.com/lissettegomez/coMethDMR/issues/6
# We met with the issue author, and we will attempt to replicate her issue in
#   order that we can identify and correct the errors in our code.

library(tidyverse)
betas <- read_delim("./inst/bug_fixes/betas.txt", delim = " ")
pheno <- read_delim("./inst/bug_fixes/pheno.txt", delim = " ")
cometh <- readRDS("./inst/bug_fixes/coMeth_16.rds")

filtered_df <- betas %>% 
	filter(CpG %in% cometh[["chr16:2961485-2962051"]]) %>% 
	as.data.frame()

filtered_df <- 
	filtered_df %>%
	column_to_rownames(var = "CpG")
filtered_df[,1:5]

# Vector of covariates
fixed_cov <- c(
	"age", "BMI", "sex", "PC1", "PC2", "PC3", "PC4", "CD8T", "CD4T", "NK",
	"Bcell", "Mono", "batch"
)

# Execute our function
library(coMethDMR)
lmmTest(
	beta = filtered_df, pheno,
	contPheno_char = "LVMindex",
	covariates_char = fixed_cov,
	modelType = "randCoef", 
	arrayType = "EPIC"
)

# Analyzing region chr16:2961485-2961799. 
# 
# chrom   start     end nCpGs Estimate StdErr pValue
# chr16 2961485 2961799     6       NA     NA      1

# TO DO: dig into the inner workings of this function to find out why both the
#   estimate and standard error are NA. Even if the algorithm fails to converge,
#   we should at least get an estimate.

