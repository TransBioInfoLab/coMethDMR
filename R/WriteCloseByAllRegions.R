

#' Extract clusters of CpG probes located closely
#'
#' @param fileName Name of the RDS file where the output genomic regions will be saved.
#' @param regions GRanges of input genomic regions
#' @param genome Human genome of reference hg19 or hg38
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param maxGap an integer, genomic locations within maxGap from each other
#'    are placed into the same cluster
#' @param minCpGs an integer, minimum number of CpGs for each resulting region
#' @param ... Dots for internal arguments. Currently unused.
#'
#' @return a file with the genomic regions containing CpGs located closely within each
#'  inputting pre-defined genomic region
#'
#' @details For \code{maxGap} = 200 and \code{minCpGs} = 3, we already calculated
#'    the clusters of CpGs. They are saved in folder \code{/inst/extdata/}.
#'
#' @export
#'
#' @importFrom dplyr group_by filter n
#' @importFrom rlang .data
#' @importFrom GenomicRanges findOverlaps
#' @examples
#'
#' regions <- GenomicRanges::GRanges(
#'   seqnames = c("chr4","chr6","chr16","chr16","chr22","chr19"),
#'   ranges = c(
#'     "174202697-174203520","28226203-28227482","89572934-89574634",
#'     "67232460-67234167","38244199-38245362","39402823-39403373"
#'   )
#' )
#'
#' \donttest{
#'   WriteCloseByAllRegions(
#'     regions = regions, arrayType = "EPIC", maxGap = 50,
#'     minCpGs = 3, fileName = "closeByRegions.rds"
#'   )
#' }
#'
WriteCloseByAllRegions <- function(
  fileName,
  regions,
  genome = c("hg19","hg38"),
  arrayType = c("450k","EPIC"),
  maxGap = 200,
  minCpGs = 3,
  ...
){

  if(maxGap == 200 & minCpGs == 3){

    warning(
      paste("A file of close by CpGs for maxGap = 200 and minCpGs = 3
            for genic and intergenic regions already exist at /inst/extdata/")
    )

  } else {

    arrayType <- match.arg(arrayType)
    genome <- match.arg(genome)

    # Available manifest files are
    # "EPIC.hg19.manifest"  "EPIC.hg38.manifest"
    # "HM27.hg19.manifest"  "HM27.hg38.manifest"
    # "HM450.hg19.manifest" "HM450.hg38.manifest"
    manifest <- paste(ifelse(arrayType == "450k","HM450","EPIC"),
                      genome, "manifest", sep = ".")
    CpGlocations.gr <- sesameDataGet(manifest)

    ### Convert input from GRanges to list of vectors of CpGs ###
    hits <- findOverlaps(regions,CpGlocations.gr) %>% as.data.frame
    hits <- hits %>%
      dplyr::group_by(.data$queryHits) %>%
      dplyr::filter(n() >= 3)
    hits$probes <- names(CpGlocations.gr)[hits$subjectHits]
    region3CpGs_ls <- split(hits$probes, hits$queryHits)

    ### Find close by clusters in all the regions ###
    ### 45.92571 secs for 1000 regions
    closeByRegions_ls <- lapply(
      region3CpGs_ls, CloseBySingleRegion, genome, arrayType, maxGap, minCpGs
    )

    ### Remove 'NULL' elements of the list ###
    closeByRegionsNoNull_ls <- unlist(closeByRegions_ls, recursive = FALSE)

    ### Order CpGs in each cluster by location to name the cluster ###
    ### 8.202557 secs for 162 regions, after unlisting 1000 regions from previous step
    closeByRegionsOrderedDF_ls <- lapply(
      closeByRegionsNoNull_ls,
      FUN = OrderCpGsByLocation,
      genome,
      arrayType,
      output = "dataframe"
    )

    ### Name each cluster with genomic region ###
    closeByRegionsNames_ls <- lapply(
      closeByRegionsOrderedDF_ls, FUN = NameRegion
    )

    ### Order CpGs in each cluster by location ###
    closeByRegionsOrdered_ls <- lapply(
      closeByRegionsOrderedDF_ls, function(x) x[,"cpg"]
    )

    names(closeByRegionsOrdered_ls) <- closeByRegionsNames_ls


    ### Select and return output ###
    message(
      paste("Writing to file", fileName)
    )

    saveRDS(closeByRegionsOrdered_ls, fileName)

  }

}
