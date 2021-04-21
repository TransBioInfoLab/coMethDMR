#' Split CpG dataframe by Subregion
#'
#' @description Split a dataframe of CpGs and comethylated subregions to a list
#'    of CpGs in each subregion
#'
#' @param CpGsSubregions_df data frame with CpG and subregion number
#' @param genome Human genome of reference hg19 or hg38
#' @param arrayType Type of array, 450k or EPIC
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
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
SplitCpGDFbyRegion <- function(CpGsSubregions_df,
                               genome = c("hg19","hg38"),
                               arrayType = c("450k", "EPIC"),
                               returnAllCpGs = TRUE){

  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  ### Output 'NULL' if thre is not a contiguous comethylated region ###
  # Function 'all' checks all the elements of the vector
  if (returnAllCpGs == FALSE & all(CpGsSubregions_df$Subregion == 0)){

    NULL

  } else {

    ### If returnAllCpGs == TRUE, output all the CpGs in the region ###

    ### Split CpGs-Subregions dataframe to list ###
    subRegion_ls <- split(
      CpGsSubregions_df$ProbeID, CpGsSubregions_df$Subregion
    )

    ### Output dataframes with annotation for each subregions ###
    subRegionAnnotationDF_ls <- lapply(
      subRegion_ls, OrderCpGsByLocation, genome, arrayType, "dataframe"
    )

    ### Name the comethylated subregions ###
    names(subRegion_ls) <- lapply(subRegionAnnotationDF_ls, NameRegion)

    subRegion_ls

  }


}
