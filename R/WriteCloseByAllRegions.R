#' Extract clusters of CpG probes located closely
#'
#' @param fileName Name of the RDS file where the output genomic regions will be
#'   saved.
#' @param regions GRanges of input genomic regions
#' @param genome Human genome of reference: hg19 or hg38
#' @param arrayType Type of array: "450k" or "EPIC"
#' @param manifest_gr A GRanges object with the genome manifest (as returned by
#'   \code{\link[ExperimentHub]{ExperimentHub}} or by
#'   \code{\link{ImportSesameData}}). This function by default ignores this
#'   argument in favour of the \code{genome} and \code{arrayType} arguments.
#' @param ignoreStrand Whether strand can be ignored, default is TRUE
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for each resulting region
#' @param ... Dots for internal arguments. Currently unused.
#'
#' @return Nothing. Instead, file with the genomic regions containing CpGs
#'   located closely within each inputting pre-defined genomic region will be 
#'   written to the disk
#'
#' @details For \code{maxGap} = 200 and \code{minCpGs} = 3, we have already
#'   calculated the clusters of CpGs. They are saved in folder
#'   \code{/inst/extdata/}.
#'
#' @export
#'
#' @importFrom dplyr group_by filter n
#' @importFrom rlang .data
#' @importFrom GenomicRanges findOverlaps
#' @importFrom methods is
#' @examples
#'
#' regions <- GenomicRanges::GRanges(
#'   seqnames = c("chr4", "chr6", "chr16", "chr16", "chr22", "chr19"),
#'   ranges = c(
#'     "174202697-174203520", "28226203-28227482", "89572934-89574634",
#'     "67232460-67234167", "38244199-38245362", "39402823-39403373"
#'   )
#' )
#'
#' # Uncomment out the example code below:
#' # WriteCloseByAllRegions(
#' #   regions = regions,
#' #   arrayType = "EPIC",
#' #   maxGap = 50,
#' #   minCpGs = 3,
#' #   fileName = "closeByRegions.rds"
#' # )
#'
#'
WriteCloseByAllRegions <- function(
  fileName,
  regions,
  genome = c("hg19","hg38"),
  arrayType = c("450k","EPIC"),
  manifest_gr = NULL,
  ignoreStrand = TRUE,
  maxGap = 200,
  minCpGs = 3,
  ...
){

  ###  Check Inputs  ###
  if (maxGap == 200 & minCpGs == 3) {
    stop(
      "A file of close by CpGs for maxGap = 200 and minCpGs = 3 for genic and
      intergenic regions already exist at /inst/extdata/",
      call. = FALSE
    )
  }
  
  if (!is(regions, "GRanges")) {
    stop(
      "The object regions must be a GRanges object from package GenomicRanges.",
      call. = FALSE
    )
  }
  
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
  
  
  ### Convert input from GRanges to list of vectors of CpGs ###
  hits <- as.data.frame( findOverlaps(regions, CpGlocations.gr) )
  # NOTE 2021-10-22: we need to re-write this line to remove the dplyr
  #   dependency. This is the only place in the whole package where we import
  #   dplyr or plyr functions.
  hits <- dplyr::filter(
    dplyr::group_by(hits, .data$queryHits),
    dplyr::n() >= 3
  )
  hits$probes <- names(CpGlocations.gr)[hits$subjectHits]
  region3CpGs_ls <- split(hits$probes, hits$queryHits)
  
  ### Find close by clusters in all the regions ###
  ### 45.92571 secs for 1000 regions
  closeByRegions_ls <- lapply(
    X = region3CpGs_ls,
    FUN = CloseBySingleRegion,
    genome,
    arrayType,
    maxGap,
    minCpGs
  )
  
  ### Remove 'NULL' elements of the list ###
  closeByRegionsNoNull_ls <- unlist(closeByRegions_ls, recursive = FALSE)
  
  ### Order CpGs in each cluster by location to name the cluster ###
  ### 8.202557 secs for 162 regions, after unlisting 1000 regions from previous step
  closeByRegionsOrderedDF_ls <- lapply(
    X = closeByRegionsNoNull_ls,
    FUN = OrderCpGsByLocation,
    genome,
    arrayType,
    manifest_gr,
    ignoreStrand,
    output = "dataframe"
  )
  
  ### Name each cluster with genomic region ###
  closeByRegionsNames_ls <- lapply(
    closeByRegionsOrderedDF_ls, NameRegion
  )
  
  ### Order CpGs in each cluster by location ###
  closeByRegionsOrdered_ls <- lapply(
    closeByRegionsOrderedDF_ls, function(x) x[ , "cpg"]
  )
  
  names(closeByRegionsOrdered_ls) <- closeByRegionsNames_ls
  
  
  ### Select and return output ###
  message("Writing to file ", fileName)
  saveRDS(closeByRegionsOrdered_ls, fileName)

}
