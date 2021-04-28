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
library(lmerTest)
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


######  Our Function Definitions  #############################################

# lmmTest
lmmTest <- function(betaOne_df, pheno_df,
										contPheno_char, covariates_char,
										modelType = c("randCoef", "simple"),
										genome = c("hg19","hg38"),
										arrayType = c("450k","EPIC"),
										outLogFile = NULL){
	# browser()
	
	modelType <- match.arg(modelType)
	arrayType <- match.arg(arrayType)
	genome <- match.arg(genome)
	
	### Transpose betaOne_df from wide to long ###
	betaOne_df$ProbeID <- row.names(betaOne_df)
	betaOneTransp_df <- reshape(
		betaOne_df,
		varying = colnames(betaOne_df)[-ncol(betaOne_df)],
		v.names = "beta",
		direction = "long",
		times = colnames(betaOne_df)[-ncol(betaOne_df)],
		timevar = "Sample"
	)
	
	### Calculate M values ###
	betaOneTransp_df$Mvalue <- log2(
		betaOneTransp_df$beta / (1 - betaOneTransp_df$beta)
	)
	
	### Merge transposed beta matrix with phenotype ###
	betaOnePheno_df <- merge(betaOneTransp_df, pheno_df, by = "Sample")
	
	# regionNames
	regionName <- NameRegion(
		OrderCpGsByLocation(
			betaOne_df$ProbeID, genome, arrayType, output = "dataframe"
		)
	)
	
	modelFormula_char <- .MakeLmmFormula(contPheno_char, covariates_char, modelType)
	
	if(!is.null(outLogFile)){
		cat(paste0("Analyzing region ", regionName, ". \n"))
	} else {
		message(paste0("Analyzing region ", regionName, ". \n"))
	}
	
	browser()
	
	f <- tryCatch({
		suppressMessages(
			lmer(as.formula(modelFormula_char), betaOnePheno_df)
		)
	}, error = function(e){NULL})
	
	
	if(is.null(f)){
		
		ps_df <- data.frame(
			Estimate = NA_real_,
			StdErr = NA_real_,
			pValue = 1
		)
		
	} else {
		
		ps_mat <- coef(summary(f))[contPheno_char, c(1, 2, 4), drop = FALSE]
		ps_df <- as.data.frame(ps_mat)
		colnames(ps_df) <- c("Estimate", "StdErr", "Stat")
		rownames(ps_df) <- NULL
		
		# If the optimization routine converged, calculate the p-value. See:
		#   https://rdrr.io/cran/lme4/man/convergence.html
		conv_ls <- f@optinfo$conv
		if(conv_ls$opt != 0){
			ps_df$pValue <- 1
		} else if(!is.null(conv_ls$lme4$messages)) {
			
			if(any(grepl("failed to converge", conv_ls$lme4$messages))){
				ps_df$pValue <- 1
			} else {
				ps_df$pValue <- 2 * (1 - pnorm(abs(ps_df$Stat)))
			}
			
		} else {
			ps_df$pValue <- 2 * (1 - pnorm(abs(ps_df$Stat)))
		}
		
	}
	
	### split regionName into chrom, start, end
	chrom <- sub(":.*",    "", regionName)
	range <- sub("c.*:",   "", regionName)
	start <- sub("-\\d*",  "", range)
	end   <- sub("\\d*.-", "", range)
	
	### Return results ###
	nCpGs <- nrow(betaOne_df)
	result <- cbind (
		chrom, start, end, nCpGs,
		ps_df,
		stringsAsFactors = FALSE
	)
	
	result
	
	
}




.MakeLmmFormula <- function(contPheno_char, covariates_char = NULL,
														modelType = c("randCoef", "simple")){
	
	modelType <- match.arg(modelType)
	
	baseMod_char <- "Mvalue ~ (1|Sample)"
	
	randomCoef_char <- paste0("(",contPheno_char, "|ProbeID)")
	
	if (!is.null(covariates_char)){
		cov_char <- paste(covariates_char, collapse = " + ")
	}
	
	######
	if(modelType == "randCoef"){
		
		ifelse(
			is.null(covariates_char),
			rcMod_char <- paste(
				baseMod_char, randomCoef_char, contPheno_char, sep = " + "
			),
			rcMod_char <- paste(
				baseMod_char, randomCoef_char, contPheno_char, cov_char, sep = " + "
			)
		)
		
		
	} else {
		
		ifelse(
			is.null(covariates_char),
			rcMod_char <- paste(baseMod_char, contPheno_char, sep = " + "),
			rcMod_char <- paste(baseMod_char, contPheno_char, cov_char, sep = " + ")
		)
		
	}
	
}


######  Try Our Function on the Test Data  ####################################

lmmTest(
	betaOne_df = filtered_df, pheno_df = pheno,
	contPheno_char = "LVMindex",
	covariates_char = fixed_cov,
	modelType = "randCoef", 
	arrayType = "EPIC"
)

# We have discovered that the original data does not have matching sample IDs.
#   A sample ID from the pheno data is "201105980120_R02C01", while a sample ID
#   from the betas data is "X200970160051_R01C01". R has set a trap for us: it
#   does not allow a base data.frame object to have column names which begin
#   with a number. Let's fix the column names, then attempt to run our code
#   again (really, this will mean pre-pending an "X" to the sample IDs in the
#   pheno data set).
pheno$Sample <- paste0("X", pheno$Sample)

# Try again:
lmmTest(
	betaOne_df = filtered_df, pheno_df = pheno,
	contPheno_char = "LVMindex",
	covariates_char = fixed_cov,
	modelType = "randCoef", 
	arrayType = "EPIC"
)

# Now we get that the object f (the output of tryCatch(suppressMessages(lmer())))
#   is an "<Object with null pointer>". Apparently, this is a "red herring" and
#   is because RStudio can't figure out what to do with an S4 object that 
#   contains an empty pointer. This is not our problem. The object "f" can still
#   be interrogated via normal class(), typeof(), str(), is.null() calls and
#   others.

# Therfore, we get these results:
# chrom   start     end nCpGs     Estimate       StdErr      Stat    pValue
# chr16 2961485 2961799     6 -0.001313689 0.0007835092 -1.676673 0.0936064
#
# It appears that the issue is sample IDs being incorrectly modified by R when
#   they are stored in the column names of a data frame. We will need to 
#   1) confirm this is the issue (re-run on a second computer independently)
#   Re-run on a second computer independently, we got the same results:
#   chrom   start     end nCpGs     Estimate       StdErr      Stat    pValue
#   1 chr16 2961485 2961799     6 -0.001313689 0.0007835092 -1.676673 0.0936064
#   2) respond with a simple fix on the GitHub issue #6
#   3) write an error in lmmTestAllRegions (and maybe even in lmmTest) to catch
#      if the sample IDs from the betas and pheno data sets do not match (that
#      is, they contain null intersection). Perhaps we can add a message if the
#      number of samples in pheno does not match the number of samples in beta.
