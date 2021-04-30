#' Extract clusters of CpGs located closely in a genomic region.
#'
#' @param CpGs_char a list of CpG IDs
#' @param genome Human genome of reference hg19 or hg38
#' @param arrayType Type of array, 450k or EPIC
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for the resulting CpG
#'    cluster
#'
#' @return a list, each item in the list is a character vector of CpG IDs
#'    located closely (i.e. in the same cluster)
#'
#' @details Note that this function depends only on CpG locations, and not on
#'    any methylation data. The algorithm is based on the
#'    \code{\link[bumphunter]{clusterMaker}} function in the \code{bumphunter}
#'    R package. Each cluster is essentially a group of CpG locations such that
#'    two consecutive locations in the cluster are separated by less than
#'    \code{maxGap}.
#'
#' @importFrom bumphunter clusterMaker
#'
#' @export
#'
#' @examples
#'
#'    CpGs_char <- c(
#'      "cg02505293", "cg03618257", "cg04421269", "cg17885402", "cg19890033",
#'      "cg20566587", "cg27505880"
#'    )
#'
#'    cluster_ls <- CloseBySingleRegion(
#'      CpGs_char,
#'      genome = "hg19",
#'      arrayType = "450k",
#'      maxGap = 100,
#'      minCpGs = 3
#'    )
#'
CloseBySingleRegion <- function(
  CpGs_char,
  genome = c("hg19","hg38"),
  arrayType = c("450k","EPIC"),
  maxGap = 200,
  minCpGs = 3
){

  CpGsOrdered_df <- OrderCpGsByLocation(
    CpGs_char, genome, arrayType, output = "dataframe"
  )

  ### Find close by clusters ###
  chr <- CpGsOrdered_df$chr
  pos <- CpGsOrdered_df$pos
  CpGsOrdered_df$cluster <- clusterMaker(chr, pos, maxGap = maxGap)

  ### Create list of vectors of CpGs in each cluster ###
  CpGsRegion_ls <- split(CpGsOrdered_df$cpg, CpGsOrdered_df$cluster)

  ### Filter for clusters with number of CpGs >= minCpGs ###
  CpGsRegionMinCpGs_ls <- CpGsRegion_ls[lengths(CpGsRegion_ls) >= minCpGs]

  if (length(CpGsRegionMinCpGs_ls) > 0 ){
    CpGsRegionMinCpGs_ls
  } else {
    NULL
  }

}
