#' Linear Mixed Models by Region
#' 
#' @description Fit mixed model to test association between a continuous
#'   phenotype and methylation values in a list of genomic regions
#'
#' @param betas data frame or matrix of beta values for all genomic regions,
#'   with row names = CpG IDs and column names = sample IDs. This is often the
#'   genome-wide array data.
#' @param region_ls a list of genomic regions; each item is a vector of CpG IDs
#'   within a genomic region. The co-methylated regions can be obtained by
#'   function \code{\link{CoMethAllRegions}}.
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'   \code{Sample} indicating sample IDs.
#' @param contPheno_char character string of the main effect (a continuous
#'   phenotype) to be tested for association with methylation values in each
#'   region
#' @param covariates_char character vector for names of the covariate variables
#' @param modelType type of mixed model; can be \code{randCoef} for random
#'   coefficient mixed model or \code{simple} for simple linear mixed model.
#' @param genome Human genome of reference: hg19 or hg38
#' @param arrayType Type of array: "450k" or "EPIC"
#' @param ignoreStrand Whether strand can be ignored, default is TRUE
#' @param outFile output .csv file with the results for the mixed model analysis
#' @param outLogFile log file for mixed models analysis messages
#' @param nCores_int Number of computing cores to be used when executing code
#'   in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'   See \code{\link{CreateParallelWorkers}} for more information.
#'
#' @return If \code{outFile} is \code{NULL}, this function returns a data frame
#'   as described below. If \code{outFile} is specified, this function writes a
#'   data frame in .csv format with the following information to the disk:
#'   location of the genomic region (\code{chrom, start, end}), number of CpGs
#'   (\code{nCpGs}), \code{Estimate}, Standard error (\code{StdErr}) of the test
#'   statistic, p-value and False Discovery Rate (FDR) for association between
#'   methylation values in each genomic region with phenotype (\code{pValue}).
#'   
#'   If \code{outLogFile} is specified, this function also writes a .txt file
#'   that includes messages for mixed model fitting to the disk.
#'
#' @details This function implements a mixed model to test association between
#'   methylation M values in a genomic region with a continuous phenotype. In
#'   our simulation studies, we found both models shown below are conservative,
#'   so p-values are estimated from normal distributions instead of Student's 
#'   \emph{t} distributions.
#'   
#'   When \code{modelType = "randCoef"}, the model is:
#'
#'   \code{M ~ contPheno_char + covariates_char + (1|Sample) + (contPheno_char|CpG)}.
#'   
#'   The last term specifies random intercept and slope for each CpG. When
#'   \code{modelType = "simple"}, the model is
#'
#'   \code{M ~ contPheno_char + covariates_char + (1|Sample)}.
#'
#'   For the results of mixed models, note that if the mixed model failed to
#'   converge, p-value for mixed model is set to 1. Also, if the mixed model is
#'   singular, at least one of the estimated variance components for intercepts
#'   or slopes random effects is 0, because there isn't enough variability in
#'   the data to estimate the random effects. In this case, the mixed model
#'   reduces to a fixed effects model. The p-values for these regions are still
#'   valid.
#'
#' @export
#'
#' @importFrom BiocParallel bplapply
#' @importFrom stats coef lm as.formula reshape p.adjust
#' @importFrom utils write.csv
#' @importFrom methods is
#'
#' @examples
#'    data(betasChr22_df)
#'    data(pheno_df)
#'
#'    CpGisland_ls <- readRDS(
#'      system.file(
#'        "extdata",
#'        "CpGislandsChr22_ex.rds",
#'        package = 'coMethDMR',
#'        mustWork = TRUE
#'      )
#'    )
#'
#'    coMeth_ls <- CoMethAllRegions(
#'      dnam = betasChr22_df,
#'      betaToM = TRUE,
#'      CpGs_ls = CpGisland_ls,
#'      arrayType = "450k",
#'      rDropThresh_num = 0.4,
#'      returnAllCpGs = FALSE
#'    )
#'
#'
#'    results_df <- lmmTestAllRegions(
#'      betas = betasChr22_df,
#'      region_ls = coMeth_ls,
#'      pheno_df = pheno_df,
#'      contPheno_char = "stage",
#'      covariates_char = "age.brain",
#'      modelType = "randCoef",
#'      arrayType = "450k",
#'      ignoreStrand = TRUE,
#'      # generates a log file in the current directory
#'      # outLogFile = paste0("lmmLog_", Sys.Date(), ".txt")
#'    )
#'
#'
lmmTestAllRegions <- function(
  betas,
  region_ls,
  pheno_df,
  contPheno_char, covariates_char,
  modelType = c("randCoef", "simple"),
  genome = c("hg19", "hg38"),
  arrayType = c("450k", "EPIC"),
  ignoreStrand = TRUE,
  outFile = NULL,
  outLogFile = NULL,
  nCores_int = 1L,
  ...){
  # browser()
  
  
  ###  Inputs  ###
  modelType <- match.arg(modelType)
  
  # Available manifest files are: 
  #   "EPIC.hg19.manifest"  "EPIC.hg38.manifest"
  #   "HM450.hg19.manifest" "HM450.hg38.manifest"
  genome <- match.arg(genome)
  arrayType <- match.arg(arrayType)
  manifest <- paste(
    switch(arrayType, "450k" = "HM450", "EPIC" = "EPIC"),
    genome, "manifest",
    sep = "."
  )
  CpGlocations.gr <- ImportSesameData(manifest)  

  if (is(betas, "matrix")){
    beta_df <- as.data.frame(betas)
  } else if (is(betas, "data.frame")) {
    beta_df <- betas
  } else {
    stop(
      "The betas object must be either a matrix or a dataframe.",
      call. = FALSE
    )
  }

  CpGnames <- rownames(beta_df)

  
  ###  Setup Log File  ###
  warnLvl <- options()$warn
  options(warn = 1)
  
  writeLog_logi <- !is.null(outLogFile)
  if(writeLog_logi){

    message("Diagnostics for mixed model fittings are in file ", outLogFile)
    # for why we need the direct file creation, two sinks opened, and two sinks
    #   closed, see the following two Stack Overflow questions:
    # 48173020/r-function-sink-isnt-redirecting-messages-or-warnings-to-a-file
    # 25948774/how-to-capture-warnings-with-the-console-output
    log_con <- file(outLogFile, open = "wt")
    sink(file = log_con)
    sink(file = log_con, type = "message")
    cat("Fitting linear mixed model to all genomic regions... \n")
    cat(paste0("Computation started at ", Sys.time(), ". \n \n"))

  }


  ###  Run mixed model for all the contiguous comethylated regions  ###
  # Split Data by Region
  coMethBetaDF_ls <- lapply(
    region_ls,
    function(x) beta_df[which(CpGnames %in% x), ]
  )
  
  # Apply
  cluster <- CreateParallelWorkers(nCores_int, ...)

  results_ls <- bplapply(
    coMethBetaDF_ls,
    FUN = lmmTest,
    BPPARAM = cluster,
    pheno_df,
    contPheno_char,
    covariates_char,
    modelType,
    genome,
    arrayType,
    manifest_gr = CpGlocations.gr,
    ignoreStrand,
    outLogFile
  )

  
  ###  Close Log File  ###
  options(warn = warnLvl)
  if(writeLog_logi){

    cat("\n")
    cat(paste0("Computation completed at ", Sys.time(), ". \n"))
    cat("Note: \n")
    cat("(1) When mixed model failed to converge, p-value for mixed model is set to 1. \n")
    cat("(2) When mixed model is singular, at least one of the estimated variance
    components for intercepts or slopes random effects is 0, because there isn't
    enough variability in data to estimate the random effects. In this case, the
    mixed model reduces to a fixed effects model. The p-values for these regions
    are still valid.\n")

    sink(type = "message")
    sink()
    close(log_con)
    # for why we need two closing sinks, see the comments above; however, even
    #   with two sink calls to close the connection, the connection is still
    #   open (it shows up with showConnections()). Thus, we also use the close()
    #   function directly.

  } else {
    message(
      "For future calls to this function, perhaps specify a log file.
      Set the file name of the log file with the outLogFile argument."
    )
  }


  ###  Output results  ###
  if (length(results_ls) > 0){

    outDF <- do.call(rbind, results_ls)
    outDF$FDR <- p.adjust(outDF$pValue, method = "fdr")
    row.names(outDF) <- NULL

  }

  if (is.null(outFile)){
    outDF
  } else {

    message("writing results to ", outFile,".csv")
    write.csv(outDF, paste0(outFile,".csv"), quote = FALSE, row.names = FALSE)

  }

}
