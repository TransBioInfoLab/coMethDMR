#' Convert genomic regions in a data frame to GRanges format
#'
#' @param regionName_char a character vector of regions, in this format: "chrxx:xxxxxx-xxxxxx"
#'
#' @return genomic regions in GRanges format
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @examples
#'  regions = c("chr22:19709548-19709755", "chr2:241721922-241722113")
#'  RegionsToRanges (regions)

RegionsToRanges <- function(regionName_char){

  chr <- sub(":.*",  "",  regionName_char)

  range <- sub ("c.*:", "",  regionName_char )

  start <- sub ("-\\d*", "", range)

  end <- sub ("\\d*.-", "", range)

  return (GRanges ( seqnames = as.factor(chr),
                    ranges = IRanges(as.numeric(start), as.numeric(end))
                  )
           )

}
