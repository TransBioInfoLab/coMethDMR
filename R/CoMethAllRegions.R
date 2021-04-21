#' Extract contiguous co-methylated genomic regions from a list of
#'   pre-defined genomic regions
#'
#' @param dnam matrix (or data frame) of beta values, with row names = CpG IDs,
#'    column names = sample IDs. This is typically genome-wide methylation beta
#'    values.
#' @param betaToM indicates if converting methylation beta values to mvalues
#' @param method method for computing correlation,
#' can be "spearman" or "pearson"
#'
#' @param rDropThresh_num threshold for min correlation between a cpg with sum
#'    of the rest of the CpGs
#' @param minCpGs minimum number of CpGs to be considered a "region".
#'    Only regions with more than \code{minCpGs} will be returned.
#'
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @param CpGs_ls list where each item is a character vector of CpGs IDs.
#'    This should be CpG probes located closely on the array.
#'
#' @param file an RDS file with clusters of CpG locations (i.e. CpGs
#'    located closely to each other on the genome). This file can be generated
#'    by the \code{\link{WriteCloseByAllRegions}} function.
#'
#' @param returnAllCpGs When there is not a contiguous comethylated region in
#'    the inputting pre-defined region, \code{returnAllCpGs = 1} indicates
#'    outputting all the CpGs in the input regions, while
#'    \code{returnAllCpGs = 0} indicates not returning any CpG.
#'
#' @param output a character vector of CpGs or a dataframe of CpGs along with
#'    rDrop info
#'
#' @param nCores_int Number of computing cores to be used when executing code
#'    in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'    See \code{\link{CreateParallelWorkers}} for more information.
#'
#' @return  When \code{output = "dataframe"} is selected,
#' returns a list of data frames, each with \code{CpG}
#' (CpG name), \code{Chr} (chromosome number), \code{MAPINFO} (genomic
#' position), \code{r_drop} (correlation between the CpG with rest of the
#' CpGs), \code{keep} (indicator for co-methylated CpG),
#' \code{keep_contiguous} (index for contiguous comethylated subregions).
#'
#' When \code{output = "CpGs"} is selected, returns a list,
#' each item is a list of CpGs in the contiguous co-methylated subregion.
#'
#' @details There are two ways to input genomic regions
#' for this function: (1) use \code{CpGs_ls} argument
#' (2) use \code{file} argument
#'
#' examples of these files are in /inst/extdata/ folder of the package.
#'
#' @export
#'
#' @importFrom BiocParallel bplapply
#'
#' @examples
#'    data(betasChr22_df)
#'
#'
#'    CpGisland_ls <- readRDS(
#'      system.file(
#'        "extdata",
#'        "CpGislandsChr22_ex.RDS",
#'        package = 'coMethDMR',
#'        mustWork = TRUE
#'      )
#'    )
#'
#'    coMeth_ls <- CoMethAllRegions (
#'      dnam = betasChr22_df,
#'      betaToM = TRUE,
#'      method = "pearson",
#'      CpGs_ls = CpGisland_ls,
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE,
#'      output = "CpGs"
#'    )
#'
#'
CoMethAllRegions <- function(dnam,
                             betaToM = FALSE,
                             method = c("pearson", "spearman"),
                             rDropThresh_num = 0.4,
                             minCpGs = 3,
                             arrayType = c("450k","EPIC"),
                             CpGs_ls,
                             file = NULL,
                             returnAllCpGs = FALSE,
                             output = c("CpGs", "dataframe"),
                             nCores_int = 1L,
                             ...){
  # browser()

  method <- match.arg(method)
  arrayType <- match.arg(arrayType)
  output <- match.arg(output)

  ### Read file of close by CpGs ###
  if(!is.null(CpGs_ls)) {

    closeByGenomicRegion_ls <- CpGs_ls

  } else {

    closeByGenomicRegion_ls <- readRDS(file)

  }

  cluster <- CreateParallelWorkers(nCores_int, ...)

  coMethCpGsAllREgions_ls <- bplapply(
    unname(closeByGenomicRegion_ls),
    FUN = CoMethSingleRegion,
    BPPARAM = cluster,
    dnam = dnam,
    betaToM = betaToM,
    rDropThresh_num = rDropThresh_num,
    minCpGs = minCpGs,
    method = method,
    arrayType = arrayType,
    returnAllCpGs = returnAllCpGs
  )

  coMethCpGsAllREgions_ls <- unique(coMethCpGsAllREgions_ls)


  ### return output ###
  # Return list of contiguous comethylated CpGs by Regions
  if(output == "CpGs"){

    unlist(
      lapply(coMethCpGsAllREgions_ls, `[[`, 2),
      recursive = FALSE
    )

  } else {

    out_ContigRegions <- lapply(coMethCpGsAllREgions_ls, `[[`, 1)
    nullRegions_lgl <- vapply(
      X = out_ContigRegions,
      FUN = is.null,
      # FUN.VALUE = logic(1); Fernanda, you **must** test your code!!
      FUN.VALUE = logical(1)
    )
    out_ContigRegions[nullRegions_lgl] <- NULL
    names(out_ContigRegions) <- unlist(lapply(out_ContigRegions, `[[`, 1, 1))

    out_ContigRegions

  }


}

