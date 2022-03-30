#' Convert genomic regions in a data frame to GRanges format
#'
#' @param regionName_char a character vector of regions in the format
#'   \code{"chrxx:xxxxxx-xxxxxx"}
#'
#' @return genomic regions in GRanges format
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @examples
#'  regions <- c("chr22:19709548-19709755", "chr2:241721922-241722113")
#'  RegionsToRanges(regions)

RegionsToRanges <- function(regionName_char){

  chr   <- sub(":.*",  "",  regionName_char)
  range <- sub(".*:", "",  regionName_char )
  start <- sub("-\\d*", "", range)
  end   <- sub("\\d*.-", "", range)

  GRanges(
    seqnames = as.factor(chr),
    ranges = IRanges(start = as.numeric(start), end = as.numeric(end))
  )

}
