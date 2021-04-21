#' Test associations of individual CpGs in a genomic region with a continuous phenotype
#'
#' @param regionName_char character string of location information for a genomic region, specified in
#' the format of "chrxx:xxxxxx-xxxxxx"
#' @param betas_df data frame of beta values with row names = CpG IDs, column names = sample IDs
#' @param pheno_df a data frame with phenotype and covariate variables, with variable "Sample" for sample IDs.
#' @param contPheno_char character string of the continuous phenotype, to be tested against methylation values
#' @param covariates_char character vector of covariate variables names
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @return a data frame with location of the genomic region (Region), CpG ID (cpg), chromosome (chr),
#' position (pos), results for testing association of methylation in individual CpGs with
#' continuous phenotype (slopeEstimate, slopePval) and annotations for the regions
#'
#' @details This function implements linear models that test association between
#' methylation values in a genomic region with a continuous phenotype. Note that methylation M values
#' are used as regression outcomes in these models. The model for each CpG is:
#'
#'   \code{methylation M value ~ contPheno_char + covariates_char}
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats reshape
#'
#' @examples
#'    data(betasChr22_df)
#'    data(pheno_df)
#'
#'    CpGsInfoOneRegion(
#'      regionName_char = "chr22:19709548-19709755",
#'      betas_df = betasChr22_df,
#'      pheno_df, contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex"),
#'      arrayType = "450k"
#'    )
#'
#'    # not adjusting for covariates
#'    CpGsInfoOneRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      betas_df = betasChr22_df,
#'      pheno_df, contPheno_char = "stage",
#'      covariates_char = NULL
#'    )
#'
CpGsInfoOneRegion <- function(regionName_char, betas_df, pheno_df,
                              contPheno_char, covariates_char,
                              arrayType = c("450k","EPIC")){

  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = {
           annotation_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
         },
         "EPIC" = {
           annotation_df = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other
         }
  )


  ### Extract individual CpGs in the region ###
  CpGsToTest_char <- GetCpGsInRegion(regionName_char, arrayType = "450k")

  ### Transpose dnam from wide to long ###
  CpGsBeta_df <- betas_df[
    which(rownames(betas_df) %in% CpGsToTest_char),
  ]

  ### Calculate M values ###
  CpGsMvalue_df <- log2(CpGsBeta_df / (1 - CpGsBeta_df))

  ### Match samples to test in pheno and beta data frames ###
  rownames(pheno_df) <- pheno_df$Sample
  samplesToTest <- intersect(colnames(CpGsMvalue_df), rownames(pheno_df))
  phenoTest_df <- pheno_df[samplesToTest, ]
  CpGsMvalueTest_df <- CpGsMvalue_df[ ,samplesToTest]
  #identical(rownames(phenoTest_df), colnames(CpGsMvalueTest_df))

  ### Run linear model for each CpG ###
  if (!is.null(covariates_char)){
    cov <- paste(covariates_char, collapse = "+")
  }

  ifelse(
    is.null(covariates_char),

    lmF <- function(Mvalue) {
      lmFormula <- as.formula(paste("Mvalue ~", contPheno_char))
      tmp = coef(summary(lm(lmFormula, data=phenoTest_df)))
      tmp[contPheno_char, c(1, 4)]
    },

    lmF <- function(Mvalue) {
      lmFormula <- as.formula(paste("Mvalue ~", contPheno_char, "+", cov))
      tmp = coef(summary(lm(lmFormula, data=phenoTest_df)))
      tmp[contPheno_char, c(1, 4)]
    }
  )

  resultAllCpGs <- data.frame(t(apply(CpGsMvalueTest_df, 1, lmF)))

  ### Return results ###
  colnames(resultAllCpGs) <- c("slopeEstimate", "slopePval")
  CpGsLocation <- OrderCpGsByLocation(
    CpGs_char = CpGsToTest_char, arrayType = arrayType, output = "dataframe"
  )
  outDF <- merge(
    CpGsLocation, resultAllCpGs,
    by.x = "cpg", by.y = "row.names", sort = FALSE
  )

  outDF$slopePval <- ifelse (
    outDF$slopePval < 0.0001,
    formatC(outDF$slopePval, format = "e", digits = 3),
    round(outDF$slopePval,4)
  )

  outDF$slopeEstimate <- round(outDF$slopeEstimate,4)

  ### Add annotations
  CpGsAnno_df <- annotation_df[CpGsToTest_char ,
                               c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")]


  outAnno_df <- merge(
    outDF,  CpGsAnno_df,
    by.x = "cpg", by.y = "row.names", sort = FALSE
  )

  outAnno_DF <- cbind(regionName_char, outAnno_df)
  colnames(outAnno_DF)[1] <- "Region"

  outAnno_DF

}


