#' Mark CpGs in contiguous and co-methylated region
#'
#' @param betaCluster_mat matrix of beta values, with rownames = sample ids,
#'    column names = CpG ids. Note that the CpGs need to be ordered by their genomic positions,
#'    this can be accomplished by the \code{OrderCpGbyLocation} function.
#'
#' @param betaToM indicates if converting to mvalues before computing correlations
#'
#' @param rDropThresh_num thershold for min correlation between a cpg with sum of the
#'    rest of the CpGs
#'
#' @param method correlation method, can be pearson or spearman
#'
#' @return A data frame with the following columns:
#'
#' \itemize{
#'   \item{\code{CpG} : }{CpG ID}
#'
#'   \item{\code{keep} : }{The CpGs with \code{keep = 1} belong to the contiguous and
#'   co-methylated region}
#'
#'   \item{\code{ind} : }{Index for the CpGs}
#'
#'   \item{\code{r_drop} : }{The correlation between each CpG with the sum of the rest of the CpGs}
#' }
#'
#' @details An outlier CpG in a genomic region will typically have low correlation with the rest of
#'  the CpGs in a genomic region. On the other hand, in a cluster of co-methylated CpGs, we expect
#'  each CpG to have high correlation with the rest of the CpGs. The \code{r.drop} statistic is used
#'  to identify these co-methylated CpGs here.
#'
#' @export
#'
#' @examples
#'    data(betaMatrix_ex1)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex1, betaToM = FALSE, method = "pearson")
#'
#'    data(betaMatrix_ex2)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex2, method = "pearson")
#'
#'    data(betaMatrix_ex3)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex3, method = "pearson")
#'
#'    data(betaMatrix_ex4)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex4, rDropThresh_num = 0.6, method = "pearson")
#'


MarkComethylatedCpGs <- function (betaCluster_mat,
                                  betaToM = TRUE,
                                  rDropThresh_num = 0.4,
                                  method = c("pearson", "spearman")) {


  ### Check that betaToM == "TRUE" only if betaCluster_mat has beta values ###

  if((min(betaCluster_mat,na.rm = TRUE) < 0 | max(betaCluster_mat > 1,na.rm = TRUE)) & betaToM == "TRUE") {
    message("The input methylation values are not beta values,
         if they are M values, 'betaToM' should be FALSE")
    return(NULL)
  }

  ### Calculate r_drop and Store CpGs ###

  if (betaToM == "TRUE") {

    mvalues_mat <- log2(betaCluster_mat / (1 - betaCluster_mat))
    clusterRdrop_df <- CreateRdrop(data = mvalues_mat, method = method)

  } else {

    clusterRdrop_df <- CreateRdrop(data = betaCluster_mat, method = method)

  }

  CpGs_char <- as.character(clusterRdrop_df$CpG)

  ### Drop CpGs with r.drop < threshold_r_num ###
  # drop these cpgs
  dropCpGs_char <- CpGs_char[clusterRdrop_df$r_drop < rDropThresh_num]

  ###  Create Output Data Frame  ###
  CpGs_df <- data.frame(
    CpG = clusterRdrop_df$CpG,
    keep = ifelse(CpGs_char %in% dropCpGs_char, 0, 1), ##(drop=0, keep=1)
    ind = seq_len(ncol(betaCluster_mat)),
    r_drop = clusterRdrop_df$r_drop,
    stringsAsFactors = FALSE
  )

  CpGs_df


}
