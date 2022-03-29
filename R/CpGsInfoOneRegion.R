#' Test Associations Between a Region and Phenotype
#' @description Test associations of individual CpGs in a genomic region with a
#'   continuous phenotype
#'
#' @param regionName_char character string of location information for a genomic
#'   region, specified in the format of \code{"chrxx:xxxxxx-xxxxxx"}
#' @param region_gr An object of class \code{\link[GenomicRanges]{GRanges}} with
#'   location information for one region. If this argument is NULL, then the 
#'   region in \code{regionName_char} is used.
#' @param betas_df data frame of beta values with row names = CpG IDs, column
#'   names = sample IDs
#' @param pheno_df a data frame with phenotype and covariate variables, with
#'   variable "Sample" for sample IDs.
#' @param contPheno_char character string of the continuous phenotype to be
#'   tested against methylation values
#' @param covariates_char character vector of covariate variables names
#' @param genome human genome of reference hg19 (default) or hg38
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param manifest_gr A GRanges object with the genome manifest (as returned by
#'   \code{\link[ExperimentHub]{ExperimentHub}} or by
#'   \code{\link{ImportSesameData}}). This function by default ignores this
#'   argument in favour of the \code{genome} and \code{arrayType} arguments.
#'
#' @return a data frame with location of the genomic region (Region), CpG ID
#'   (cpg), chromosome (chr), position (pos), results for testing association of
#'   methylation in individual CpGs with continuous phenotype (slopeEstimate,
#'   slopePval) and annotations for the region.
#'
#' @details This function implements linear models that test association between
#'   methylation values in a genomic region with a continuous phenotype. Note
#'   that methylation M values are used as regression outcomes in these models.
#'   The model for each CpG is:
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
#'    myRegion_gr <- RegionsToRanges("chr22:18267969-18268249")
#'
#'    CpGsInfoOneRegion(
#'      region_gr = myRegion_gr,
#'      betas_df = betasChr22_df,
#'      pheno_df = pheno_df,
#'      contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex"),
#'      arrayType = "450k"
#'    )
#'
CpGsInfoOneRegion <- function(
  regionName_char,
  region_gr = NULL,
  betas_df,
  pheno_df,
  contPheno_char,
  covariates_char = NULL,
  genome = c("hg19", "hg38"),
  arrayType = c("450k", "EPIC"),
  manifest_gr = NULL
){
  # browser()

  
  ###  Inputs  ###
  genome <- match.arg(genome)
  arrayType <- match.arg(arrayType)

  switch(
    arrayType,
    "450k" = {
      annotation_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
    },
    "EPIC" = {
      annotation_df = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other
    }
  )


  ### Extract individual CpGs in the region ###
  if (is.null(region_gr)) {
    
    CpGsToTest_char <- GetCpGsInRegion(
      regionName_char = regionName_char,
      genome = genome,
      arrayType = arrayType, 
      manifest_gr = manifest_gr
    )
    
  } else {
    
    CpGsToTest_char <- GetCpGsInRegion(
      region_gr = region_gr,
      genome = genome,
      arrayType = arrayType, 
      manifest_gr = manifest_gr
    )
    regionName_char <- as.character(region_gr)
    
  }

  
  ###  Wrangle and Tidy Data  ###
  CpGsBeta_df <- betas_df[rownames(betas_df) %in% CpGsToTest_char, ]
  CpGsMvalue_df <- log2(CpGsBeta_df / (1 - CpGsBeta_df))

  # Match samples to test in pheno and beta data frames
  rownames(pheno_df) <- pheno_df$Sample
  samplesToTest <- intersect(colnames(CpGsMvalue_df), rownames(pheno_df))
  phenoTest_df <- pheno_df[samplesToTest, ]
  CpGsMvalueTest_df <- CpGsMvalue_df[ , samplesToTest]

  
  ###  Function to run linear model for each CpG  ###
  if (is.null(covariates_char)){
    
    lmF <- function(Mvalue) {
      lmFormula <- as.formula(paste("Mvalue ~", contPheno_char))
      tmp = coef(summary(lm(lmFormula, data = phenoTest_df)))
      tmp[contPheno_char, c(1, 4)]
    }
    
  } else {
    
    cov <- paste(covariates_char, collapse = "+")
    lmF <- function(Mvalue) {
      lmFormula <- as.formula(paste("Mvalue ~", contPheno_char, "+", cov))
      tmp <- coef(summary(lm(lmFormula, data = phenoTest_df)))
      tmp[contPheno_char, c(1, 4)]
    }
    
  }

  
  ###  Run the Models  ###
  resultAllCpGs <- data.frame(
    t( apply(CpGsMvalueTest_df, 1, lmF) )
  )
  colnames(resultAllCpGs) <- c("slopeEstimate", "slopePval")
  

  ###  Wrangle results  ###
  CpGsLocation <- OrderCpGsByLocation(
    CpGs_char = CpGsToTest_char,
    genome = genome,
    arrayType = arrayType,
    manifest_gr = manifest_gr,
    output = "dataframe"
  )
  
  outDF <- merge(
    CpGsLocation, resultAllCpGs,
    by.x = "cpg", by.y = "row.names", sort = FALSE
  )

  outDF$slopePval <- ifelse (
    outDF$slopePval < 0.0001,
    formatC(outDF$slopePval, format = "e", digits = 3),
    round(outDF$slopePval, 4)
  )

  outDF$slopeEstimate <- round(outDF$slopeEstimate, 4)

  
  ###  Add Annotations and Return  ###
  CpGsAnno_df <- annotation_df[
    CpGsToTest_char,
    c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")
  ]

  outAnno_df <- merge(
    outDF, CpGsAnno_df,
    by.x = "cpg", by.y = "row.names", sort = FALSE
  )

  outAnno_DF <- cbind(regionName_char, outAnno_df)
  colnames(outAnno_DF)[1] <- "Region"

  outAnno_DF

}


