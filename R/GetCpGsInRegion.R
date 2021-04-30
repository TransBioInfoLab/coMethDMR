
#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param arrayType Type of array, 450k or EPIC
#' @param genome human genome of reference hg19 (default) or hg38
#'
#' @return vector of CpG probe IDs mapped to the genomic region
#' @export
#'
#' @importFrom tidyr separate %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#'    GetCpGsInRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      genome = "hg19",
#'      arrayType = "450k"
#'    )
GetCpGsInRegion <- function(
  regionName_char,
  genome = c("hg19","hg38"),
  arrayType = c("450k","EPIC")
){

  arrayType <- match.arg(arrayType)
  genome <- match.arg(genome)

  # Available manifest files are
  # "EPIC.hg19.manifest"  "EPIC.hg38.manifest"
  # "HM27.hg19.manifest"  "HM27.hg38.manifest"
  # "HM450.hg19.manifest" "HM450.hg38.manifest"
  manifest <- paste(
    ifelse(arrayType == "450k","HM450","EPIC"),
    genome, "manifest",
    sep = "."
  )
  CpGlocations.gr <- sesameDataGet(manifest)

  ### Split the region name in chr and positions ###
  gr <- regionName_char %>%
    as.data.frame %>%
    separate(col = ".", into = c("seqnames","start","end")) %>%
    makeGRangesFromDataFrame()

  CpGlocations.gr <- subsetByOverlaps(CpGlocations.gr, gr)

  OrderCpGsByLocation(
    names(CpGlocations.gr),
    genome = genome,
    arrayType = arrayType,
    output = "vector"
  )
}
