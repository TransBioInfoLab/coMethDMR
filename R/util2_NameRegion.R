
#' Name a region with several CpGs based on its genomic location
#'
#' @param CpGsOrdered_df dataframe with columns for Probe IDs as character (cpg),
#'    chromosome number as character (chr) and genomic location as integer (pos)
#'
#' @return genome location of the CpGs, in the format of "chrxx:xxxxxx-xxxxxx"
#' @export
#'
#' @examples
#'   CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#'   CpGsOrdered_df <- OrderCpGsByLocation(CpGs_char, arrayType=c("EPIC"), output = "dataframe")
#'   NameRegion(CpGsOrdered_df)
#'
NameRegion <- function(CpGsOrdered_df){

  ### Return region name based on genomic location ###
    paste0(
    CpGsOrdered_df$chr[1], ":",
    CpGsOrdered_df$pos[1], "-", CpGsOrdered_df$pos[nrow(CpGsOrdered_df)]
  )


}
