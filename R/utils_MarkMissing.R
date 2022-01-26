#' Return Column and Row Names of Samples and Probes under the Missingness
#'   Theshold
#'
#' @param dnaM_df A data frame of DNA methylation values. Samples are columns.
#'   Row names are probe IDs.
#' @param sampMissing_p The maximum proportion of missingness allowed in a 
#'   sample. Defaults to 50\%.
#' @param probeMissing_p The maximum proportion of missingness allowed in a
#'   probe. Defaults to 25\%.
#'
#' @return A list of four entries: 
#' \itemize{
#'   \item{\code{dropSamples}: }{the column names of samples with more than
#'     \code{sampMissing_p} percent missing values}
#'   \item{\code{keepSamples}: }{the column names of samples with less than or
#'     equal to \code{sampMissing_p} percent missing values}
#'   \item{\code{dropProbes}: }{the row names of probes with more than
#'     \code{probeMissing_p} percent missing values}
#'   \item{\code{keepProbes}: }{the row names of probes with less than or equal
#'     to \code{probeMissing_p} percent missing values}
#' }
#' 
#' @details Before calculating the missing proportion of samples, probes with
#'   missingness greater than the threshold are dropped first. 
#' 
#' @export
#'
#' @examples
#' 
#'   ###  Setup  ###
#'   values_num <- c(
#'     0.1, 0.1, 0.1, 0.1, 0.1,
#'     0.1, 0.1, 0.1, 0.1,  NA,
#'     0.1, 0.1, 0.1, 0.1,  NA,
#'     0.1, 0.1, 0.1,  NA,  NA,
#'     0.1, 0.1, 0.1,  NA,  NA,
#'     0.1, 0.1,  NA,  NA,  NA,
#'     0.1, 0.1,  NA,  NA,  NA,
#'     0.1,  NA,  NA,  NA,  NA,
#'      NA,  NA,  NA,  NA,  NA
#'   )
#'   values_mat <- matrix(values_num, nrow = 9, ncol = 5, byrow = TRUE)
#'   rownames(values_mat) <- paste0("probe_0", 1:9)
#'   colnames(values_mat) <- paste0("sample_0", 1:5)
#'   values_df <- as.data.frame(values_mat)
#'   
#'   
#'   ###  Simple Calculations  ###
#'   MarkMissing(values_df)
#'   MarkMissing(values_df, probeMissing_p = 0.5)
#'   MarkMissing(values_df, sampMissing_p = 0.25)
#'   
#'   
#'   ###  Using the Output  ###
#'   mark_ls <- MarkMissing(values_df, probeMissing_p = 0.5)
#'   valuesPurged_df <- values_df[ mark_ls$keepProbes, mark_ls$keepSamples ]
#'   valuesPurged_df
#'   
MarkMissing <- function(dnaM_df, sampMissing_p = 0.5, probeMissing_p = 0.25) {
  # browser()
  
  # Helper function for missing proportion
  propMissing <- function(x) { mean( is.na(x) ) }
  
  # Calculate "good" probes
  goodProbes_lgl  <- apply(dnaM_df, MARGIN = 1, propMissing) <= probeMissing_p
  badProbes_char  <- row.names(dnaM_df)[!goodProbes_lgl]
  goodProbes_char <- row.names(dnaM_df)[goodProbes_lgl]
  # Keep "good" probes
  dnaM_df <- dnaM_df[goodProbes_lgl, ]
  
  # Calculate "good" samples (after removing "bad" probes)
  goodSamps_lgl <- 
    vapply(dnaM_df, propMissing, FUN.VALUE = numeric(1)) <= sampMissing_p
  
  # Return
  list(
    dropSamples = colnames(dnaM_df)[!goodSamps_lgl],
    keepSamples = colnames(dnaM_df)[goodSamps_lgl],
    dropProbes  = badProbes_char,
    keepProbes  = goodProbes_char
  )
  
}
