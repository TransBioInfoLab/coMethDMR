#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one
#'   region in the format \code{"chrxx:xxxxxx-xxxxxx"}
#' @param region_gr An object of class \code{\link[GenomicRanges]{GRanges}} with
#'   location information for one region. If this argument is NULL, then the 
#'   region in \code{regionName_char} is used.
#' @param genome human genome of reference hg19 (default) or hg38
#' @param arrayType Type of array, 450k or EPIC
#' @param manifest_gr A GRanges object with the genome manifest (as returned by
#'   \code{\link[ExperimentHub]{ExperimentHub}} or by
#'   \code{\link{ImportSesameData}}). This function by default ignores this
#'   argument in favour of the \code{genome} and \code{arrayType} arguments.
#' @param ignoreStrand Whether strand can be ignored, default is TRUE
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
#'      arrayType = "450k",
#'      ignoreStrand = TRUE
#'    )
#'    
GetCpGsInRegion <- function(
  regionName_char,
  region_gr = NULL,
  genome = c("hg19", "hg38"),
  arrayType = c("450k", "EPIC"),
  manifest_gr = NULL,
  ignoreStrand = TRUE
){
  
  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)
  
  
  ###  The GRanges Object  ###
  if(!is.null(manifest_gr)) {
    CpGlocations.gr <- manifest_gr
  } else {
    
    manifest <- paste(
      switch(arrayType, "450k" = "HM450", "EPIC" = "EPIC"),
      genome, "manifest",
      sep = "."
    )
    CpGlocations.gr <- ImportSesameData(manifest)  
    
  }
  
  
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
    manifest_gr = manifest_gr,
    ignoreStrand = ignoreStrand,
    output = "vector"
  )
  
}
