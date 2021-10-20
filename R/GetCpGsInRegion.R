#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one
#'   region in the format \code{"chrxx:xxxxxx-xxxxxx"}
#' @param region_gr An object of class \code{\link[GenomicRanges]{GRanges}} with
#'   location information for one region. If this argument is NULL, then the 
#'   region in \code{regionName_char} is used.
#' @param arrayType Type of array, 450k or EPIC
#' @param genome human genome of reference hg19 (default) or hg38
#'
#' @return vector of CpG probe IDs mapped to the genomic region
#' 
#' @export
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#'    
#'    myRegion_gr <- RegionsToRanges("chr22:18267969-18268249")
#'    
#'    GetCpGsInRegion(
#'      region_gr = myRegion_gr,
#'      genome = "hg19",
#'      arrayType = "450k"
#'    )
#'    
GetCpGsInRegion <- function(
  regionName_char,
  region_gr = NULL,
  genome = c("hg19", "hg38"),
  arrayType = c("450k", "EPIC")
){

  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  # Available manifest files are
  # "EPIC.hg19.manifest"  "EPIC.hg38.manifest"
  # "HM27.hg19.manifest"  "HM27.hg38.manifest"
  # "HM450.hg19.manifest" "HM450.hg38.manifest"
  manifest <- paste(
    ifelse(arrayType == "450k", "HM450", "EPIC"),
    genome, "manifest",
    sep = "."
  )
  CpGlocations.gr <- sesameDataGet(manifest)

  ### Split the region name in chr and positions ###
  if (is.null(region_gr)) {
    gr <- RegionsToRanges(regionName_char)
  } else {
    gr <- region_gr
  }
  CpGlocations.gr <- subsetByOverlaps(CpGlocations.gr, gr)

  OrderCpGsByLocation(
    names(CpGlocations.gr),
    genome = genome,
    arrayType = arrayType,
    output = "vector"
  )
  
}
