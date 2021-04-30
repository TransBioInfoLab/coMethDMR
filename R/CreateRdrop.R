#' Computes leave-one-out correlations (rDrop) for each CpG
#'
#' @param data a dataframe with rownames = sample IDs, column names = CpG IDs.
#' @param method method for computing correlation, can be "pearson" or "spearman"
#'
#' @importFrom stats cor
#'
#' @return A data frame with the following columns:
#'
#' \itemize{
#'   \item{\code{CpG} : }{CpG ID}
#'   \item{\code{r_drop} : }{The correlation between each CpG with the sum of
#'   the rest of the CpGs}
#' }
#' @details An outlier CpG in a genomic region will typically have low correlation with the rest of
#'  the CpGs in a genomic region. On the other hand, in a cluster of co-methylated CpGs, we expect
#'  each CpG to have high correlation with the rest of the CpGs. The \code{r.drop} statistic is used
#'  to identify these co-methylated CpGs here.
#'
#'
#' @export
#'
#' @examples
#'    data(betaMatrix_ex1)
#'
#'    CreateRdrop(data = betaMatrix_ex1, method = "pearson")
#'
CreateRdrop <- function(data, method = c("pearson", "spearman")){
  # browser()

  method <- match.arg(method)
  col_data <- colnames(data)

  out_num <- vapply(
    X = seq_along(col_data),
    FUN = function(column){

      ## remove site i and then compute row mean
      data_no_i <- data[, -column, drop = FALSE]
      data_i <- data[, column]

      data_no_i_mean <- rowMeans(data_no_i)

      ## Correlate dat_i and dat_no_i_mean
      cor(data_i, data_no_i_mean, method = method)

    },
    FUN.VALUE = numeric(1)
  )

  data.frame(
    CpG = col_data,
    r_drop = out_num,
    stringsAsFactors = FALSE
  )

}


