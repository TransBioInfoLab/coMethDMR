#' Find Contiguous Co-Methylated Regions
#' 
#' @description Find contiguous comethylated regions based on output file from
#'   function \code{MarkComethylatedCpGs}
#'
#' @param CpGs_df an output dataframe from function \code{MarkComethylatedCpGs},
#'   with variables: \code{CpG, keep, ind, r_drop}. See details in documentation
#'   for \code{\link{MarkComethylatedCpGs}}.
#' @param minCpGs_int an integer indicating the minimum number of CpGs for
#'   output genomic regions
#'
#' @return A data frame with variables \code{ProbeID} and \code{Subregion}
#'   (index for each output contiguous comethylated region)
#'
#' @export
#'
#' @importFrom bumphunter getSegments
#' @importFrom utils globalVariables
#'
#' @examples
#'    data(betaMatrix_ex4)
#'
#'    CpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex4)
#'
#'    FindComethylatedRegions(CpGs_df)
#'
FindComethylatedRegions <- function(CpGs_df, minCpGs_int = 3){
  # browser()

  ### Get contiguous regions of CpGs ###
  contiguousRegion_ls <- getSegments(x = CpGs_df$keep, cutoff = 1)[["upIndex"]]
  nSegs_int <- length(contiguousRegion_ls)

  if (nSegs_int == 0){
    
    contiguousRegionsCpGs <- CpGs_df["CpG"]
    contiguousRegionsCpGs["Subregion"] <- 0
    
  } else {

    ### Select segments with number of CpGs >= minCpGs ###
    contiguous_int <- lengths(contiguousRegion_ls)
    contiguousMinCpGs_idx <- which(contiguous_int >= minCpGs_int)
    nSegsMinCpGs_int <- length(contiguousMinCpGs_idx)

    ### Create output dataframe ###
    # Keep CpGs and contiguous comethylated subregion number
    ind <- NULL
    if (nSegsMinCpGs_int == 0){
      
      contiguousRegionsCpGs <- CpGs_df["CpG"]
      contiguousRegionsCpGs["Subregion"] <- 0
      
    } else {
      
      burner_fun <- function(u){
        
        regionSet <- contiguousRegion_ls[[contiguousMinCpGs_idx[u]]]
        whichCpG <- subset(CpGs_df, ind %in% regionSet, select = "CpG")
        subRegion_idx <- rep(u, length(regionSet))
        
        data.frame(CpG = whichCpG, subregion = subRegion_idx)
        
      }
      
      inner_ls <- lapply(seq_len(nSegsMinCpGs_int), burner_fun)
      contiguousRegionsCpGs <- do.call(rbind, inner_ls)
      
    } 

  } 

  colnames(contiguousRegionsCpGs) <- c("ProbeID", "Subregion")

  contiguousRegionsCpGs


}
