#' Get Linear Model Residuals
#' 
#' @description Remove covariate effects from methylayion values by fitting
#'   probe-specific linear models
#'
#' @param dnam data frame or matrix of methylation values with row names = CpG
#'   IDs and column names = sample IDs. This is often the genome-wide array
#'   data. 
#' @param betaToM indicates if methylation beta values (ranging from [0, 1])
#'   should be converted to M values (ranging from (-Inf, Inf)). Note that if
#'   beta values are the input to \code{dnam}, then \code{betaToM} should be set
#'   to \code{TRUE}, otherwise \code{FALSE}.
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'   \code{Sample} indicating sample IDs.
#' @param covariates_char character vector for names of the covariate variables
#' @param nCores_int Number of computing cores to be used when executing code
#'   in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'   See \code{\link{CreateParallelWorkers}} for more information.
#'    
#' @details This function fits an ordinary linear model predicting methylation
#'   values for each probe from the specified covariates. This process will be
#'   useful in scenarios where methylation values in a region or at an
#'   individual probe are known \emph{a priori} to have differential methylation
#'   independent of the disease or condition of interest.
#'
#' @return output a matrix of residual values in the same dimension as
#'   \code{dnam}
#'
#' @export
#'
#' @importFrom stats na.exclude
#' @importFrom stats residuals
#' @importFrom methods is
#'
#' @examples
#'    data(betasChr22_df)
#'
#'    data(pheno_df)
#'
#'    GetResiduals(
#'      dnam = betasChr22_df[1:10, 1:10],
#'      betaToM = TRUE,
#'      pheno_df = pheno_df,
#'      covariates_char = c("age.brain", "sex", "slide")
#'    )
#'
GetResiduals <- function(
  dnam,
  betaToM = TRUE,
  pheno_df,
  covariates_char,
  nCores_int = 1L,
  ...
){
  
  # browser()

  if (is(dnam, "matrix")){
    dnam_df = as.data.frame(dnam)
  } else {
    dnam_df = dnam
  }


  if (betaToM){
    ### Compute M values
    value_df <- log2(dnam_df / (1 - dnam_df))
  } else {
    value_df <- dnam_df
  }

  ### Select samples in both value_df and pheno_df
  if(!"Sample" %in% colnames(pheno_df)) {
    message("Could not find sample column in the pheno_df.")
    return(NULL)
  }

  ## Make sure Sample is character but not factor
  pheno_df$Sample <- as.character(pheno_df$Sample)

  ## Check if samples in value_df and pheno_df are identical
  idt <- identical(colnames(value_df), pheno_df$Sample)

  if (idt) {

    value_df <- value_df
    pheno_df <- pheno_df

  } else {
    
    message("Phenotype data is not in the same order as methylation data. We will use column Sample in phenotype data to put these two files in the same order.")
    intersectSample <- intersect(colnames(value_df), pheno_df$Sample)

    ### Select samples of pheno_df based on intersect samples
    pheno_df <- pheno_df[pheno_df$Sample %in% intersectSample, ]

    ### Select samples of value_df in pheno_df
    value_df <- value_df[ , pheno_df$Sample]

  }

  ### Create the formula
  cov_char <- paste(covariates_char, collapse = " + ")
  formula_char <- paste0("val ~ ", cov_char)

  cluster <- CreateParallelWorkers(nCores_int, ...)

  resid_ls <- bplapply(
    seq_len(nrow(value_df)),
    function(row){

      val <- t(value_df[row, ])
      colnames(val) <- "val"

      dat <- cbind(val, pheno_df)
      dat$val <- as.numeric(dat$val)

      fitE <- lm(formula_char, data = dat, na.action = na.exclude)

      residuals(fitE)

    },
    BPPARAM = cluster
  )
  
  ### Take residuals
  resid_mat <- do.call(rbind, resid_ls)
  row.names(resid_mat) <- row.names(value_df)

  resid_mat

}
