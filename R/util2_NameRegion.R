#' Name a region with several CpGs based on its genomic location
#'
#' @param CpGsOrdered_df dataframe with columns for Probe IDs as character
#'   (cpg), chromosome number as character (chr), and genomic location as
#'   integer (pos)
#'
#' @return genome location of the CpGs formatted as \code{"chrxx:xxxxxx-xxxxxx"}
#' 
#' @export
#'
#' @examples
#'  # Consider four probe IDs:
#'  CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#'  
#'  # After querying these four probes against an EPIC array via the 
#'  #   OrderCpGsByLocation() function, we get the following data frame:
#'  CpGsOrdered_df <- data.frame(
#'    chr = c("chr10", "chr10", "chr10", "chr10"),
#'    pos = c(100028236L, 100028320L, 100028468L, 100028499L),
#'    cpg = c("cg20214853", "cg04677227", "cg11632906", "cg07146435"),
#'    stringsAsFactors = FALSE
#'  )
#' 
#'  # Now, we can name the region that contains these four probes:
#'  NameRegion(CpGsOrdered_df)
#' 
NameRegion <- function(CpGsOrdered_df){
  
  range_char <- c(
    CpGsOrdered_df$pos[1],
    CpGsOrdered_df$pos[nrow(CpGsOrdered_df)]
  )
  
  paste0( CpGsOrdered_df$chr[1], ":", range_char[1], "-", range_char[2] )
  
}
