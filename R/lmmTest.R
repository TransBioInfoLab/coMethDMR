#' Fit mixed model to methylation values in one genomic region
#'
#' @param betaOne_df matrix of beta values for one genomic region, with row
#'   names = CpG IDs and column names = sample IDs
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'   \code{Sample} indicating sample IDs.
#' @param contPheno_char character string of the main effect (a continuous
#'   phenotype) to be tested for association with methylation values in the
#'   region
#' @param covariates_char character vector for names of the covariate variables
#' @param modelType type of mixed model: can be \code{randCoef} for random
#'    coefficient mixed model or \code{simple} for simple linear mixed model.
#' @param genome Human genome of reference: hg19 or hg38
#' @param arrayType Type of array: "450k" or "EPIC"
#' @param manifest_gr A GRanges object with the genome manifest (as returned by
#'   \code{\link[ExperimentHub]{ExperimentHub}} or by
#'   \code{\link{ImportSesameData}}). This function by default ignores this
#'   argument in favour of the \code{genome} and \code{arrayType} arguments.
#' @param ignoreStrand Whether strand can be ignored, default is TRUE
#' @param outLogFile Name of log file for messages of mixed model analysis
#'
#' @return  A dataframe with one row for association result of one region and 
#'   the following columns: \code{Estimate}, \code{StdErr}, and \code{pvalue}
#'   showing the association of methylation values in the genomic region tested
#'   with the continuous phenotype supplied in \code{contPheno_char}
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
#' @export
#'
#' @importFrom lmerTest lmer
#' @importFrom stats coef pnorm reshape
#'
#' @examples
#'   data(betasChr22_df)
#'
#'   CpGsChr22_char <- c(
#'     "cg02953382", "cg12419862", "cg24565820", "cg04234412", "cg04824771",
#'     "cg09033563", "cg10150615", "cg18538332", "cg20007245", "cg23131131",
#'     "cg25703541"
#'   )
#'
#'   coMethCpGs <- CoMethSingleRegion(CpGsChr22_char, betasChr22_df)
#'
#'   # test only the first co-methylated region
#'   coMethBeta_df <- betasChr22_df[coMethCpGs$CpGsSubregions[[1]], ]
#'
#'   data(pheno_df)
#'
#'   res <- lmmTest(
#'     betaOne_df = coMethBeta_df,
#'     pheno_df,
#'     contPheno_char = "stage",
#'     covariates_char = c("age.brain", "sex"),
#'     modelType = "randCoef",
#'     arrayType = "450k", 
#'     ignoreStrand = TRUE
#'   )
#'

lmmTest <- function(betaOne_df, pheno_df, contPheno_char, covariates_char,
                    modelType = c("randCoef", "simple"),
                    genome = c("hg19", "hg38"),
                    arrayType = c("450k", "EPIC"),
                    manifest_gr = NULL,
                    ignoreStrand = TRUE,
                    outLogFile = NULL){
  # browser()

  ###  Inputs  ###
  modelType <- match.arg(modelType)
  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  
  ###  Wrangle  ###
  betaOne_df$ProbeID <- row.names(betaOne_df)
  betaOneTransp_df <- reshape(
    betaOne_df,
    varying = colnames(betaOne_df)[-ncol(betaOne_df)],
    v.names = "beta",
    direction = "long",
    times = colnames(betaOne_df)[-ncol(betaOne_df)],
    timevar = "Sample"
  )

  # Calculate M values
  betaOneTransp_df$Mvalue <- log2(
    betaOneTransp_df$beta / (1 - betaOneTransp_df$beta)
  )

  # Merge transposed beta matrix with phenotype
  betaOnePheno_df <- merge(betaOneTransp_df, pheno_df, by = "Sample")

  
  ###  Setup  ###
  regionName <- NameRegion(
    OrderCpGsByLocation(
      CpGs_char = betaOne_df$ProbeID,
      genome = genome,
      arrayType = arrayType,
      manifest_gr = manifest_gr,
      ignoreStrand = ignoreStrand,
      output = "dataframe"
    )
  )
  
  if(!is.null(outLogFile)){
    cat(paste0("Analyzing region ", regionName, ". \n"))
  } else {
    message("Analyzing region ", regionName, ". \n")
  }

  modelFormula_char <- .MakeLmmFormula(
    contPheno_char = contPheno_char,
    covariates_char = covariates_char,
    modelType = modelType
  )

  
  ###  Analysis  ###
  # Run
  f <- tryCatch(
    {
      suppressMessages(
        lmer(as.formula(modelFormula_char), betaOnePheno_df)
      )
    },
    error = function(e){NULL}
  )

  # Check
  if(is.null(f)){

    ps_df <- data.frame(
      Estimate = NA_real_,
      StdErr = NA_real_,
      Stat = NA_real_,
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

  
  ###  Return  ###
  # split regionName into chrom, start, end
  chrom <- sub(":.*",    "", regionName)
  range <- sub("c.*:",   "", regionName)
  start <- sub("-\\d*",  "", range)
  end   <- sub("\\d*.-", "", range)

  # results
  nCpGs <- nrow(betaOne_df)
  cbind(
    chrom, start, end, nCpGs, ps_df,
    stringsAsFactors = FALSE
  )

}



.MakeLmmFormula <- function(contPheno_char, covariates_char = NULL,
                            modelType = c("randCoef", "simple")){

  modelType <- match.arg(modelType)

  baseMod_char <- "Mvalue ~ (1|Sample)"
  randomCoef_char <- paste0("(",contPheno_char, "|ProbeID)")

  if (!is.null(covariates_char)){
    cov_char <- paste(covariates_char, collapse = " + ")
  }

  
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
