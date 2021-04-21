

#' Find contiguous comethylated regions based on output file from function \code{MarkComethylatedCpGs}
#'
#' @param CpGs_df an output dataframe from function \code{MarkComethylatedCpGs}, with variables
#'
#'  \code{CpG, keep, ind, r_drop}. See details in documentation for \code{MarkComethylatedCpGs}.
#'
#' @param minCpGs_int an integer, indicates minimum number of CpGs for output genomic regions
#'
#' @return A data frame with variables \code{ProbeID} and \code{Subregion} (index for each output
#' contiguous comethylated regions)
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

  ### Get contiguous regions of CpGs ###
  contiguousRegion_ls <- getSegments(CpGs_df$keep, cutoff = 1)
  nSegs_int <- length(contiguousRegion_ls$upIndex)

  if (nSegs_int > 0){

    ### Select segments with number of CpGs >= minCpGs ###
    contiguous_int <- lengths(contiguousRegion_ls$upIndex)
    contiguousMinCpGs_idx <- which(contiguous_int >= minCpGs_int)
    nSegsMinCpGs_int <- length(contiguousMinCpGs_idx)

    ### Create output dataframe with CpGs and contiguous comethylated subregion number ###
    #globalVariables("ind")
    ind <- NULL

       if (nSegsMinCpGs_int > 0){

         inner_ls <- lapply(seq_len(nSegsMinCpGs_int), function(u){

           data.frame(
             CpG = subset(
               CpGs_df,
               ind %in% contiguousRegion_ls$upIndex[[contiguousMinCpGs_idx[u]]],
               select = "CpG"
             ),
             subregion = rep(
               u, length(contiguousRegion_ls$upIndex[[contiguousMinCpGs_idx[u]]])
             )
           )

         })

         contiguousRegionsCpGs <- do.call(rbind, inner_ls)


       } else {
             contiguousRegionsCpGs <- cbind(
               as.data.frame(CpGs_df$CpG),
               rep(0,length(CpGs_df$CpG))
             )
         }


  } else {
      contiguousRegionsCpGs <- cbind(
        as.data.frame(CpGs_df$CpG),
        rep(0,length(CpGs_df$CpG))
      )
    }

  colnames(contiguousRegionsCpGs) <- c("ProbeID","Subregion")

  contiguousRegionsCpGs


}
