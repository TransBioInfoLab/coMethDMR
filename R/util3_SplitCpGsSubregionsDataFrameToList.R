#' Split CpG dataframe by Subregion
#'
#' @description Split a dataframe of CpGs and comethylated subregions to a list
#'   of CpGs in each subregion
#'
#' @param CpGsSubregions_df data frame with CpG and subregion number
#' @param genome Human genome of reference: hg19 or hg38
#' @param arrayType Type of array: 450k or EPIC
#' @param manifest_gr A GRanges object with the genome manifest (as returned by
#'   \code{\link[ExperimentHub]{ExperimentHub}} or by
#'   \code{\link{ImportSesameData}}). This function by default ignores this
#'   argument in favour of the \code{genome} and \code{arrayType} arguments.
#' @param returnAllCpGs indicates if outputting all the CpGs in the region when
#'   there is not a contiguous comethylated region or only the CpGs in the
#'   contiguous comethylated regions
#'
#' @return a list of comethylated subregions CpGs for a pre-defined region
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'    data(betaMatrix_ex4)
#'    CpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex4)
#'    CpGsSubregions_df <- FindComethylatedRegions(CpGs_df)
#'
#'    SplitCpGDFbyRegion(
#'      CpGsSubregions_df,
#'      genome = "hg19",
#'      arrayType = "450k"
#'    )
#'
SplitCpGDFbyRegion <- function(
  CpGsSubregions_df,
  genome = c("hg19", "hg38"),
  arrayType = c("450k", "EPIC"),
  manifest_gr = NULL,
  returnAllCpGs = TRUE
){

  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  if (returnAllCpGs == FALSE & all(CpGsSubregions_df$Subregion == 0)){
    # Output 'NULL' if there is not a contiguous comethylated region
    return(NULL)
  } 
  
  # If returnAllCpGs == TRUE, or if there is more than one subregion, output all
  #   the CpGs in those regions
  
  # Split CpGs-subregions dataframe to list
  subRegion_ls <- split(
    CpGsSubregions_df$ProbeID, CpGsSubregions_df$Subregion
  )
  
  # Output dataframes with annotation for each subregions
  subRegionAnnotationDF_ls <- lapply(
    X = subRegion_ls,
    FUN = OrderCpGsByLocation,
    genome,
    arrayType,
    manifest_gr,
    output = "dataframe"
  )
  
  ### Name the comethylated subregions ###
  names(subRegion_ls) <- lapply(subRegionAnnotationDF_ls, NameRegion)
  
  subRegion_ls

}
