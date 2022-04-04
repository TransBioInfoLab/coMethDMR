#' @title Import Illumina manifests (sesameData versions)
#' 
#' @description Load either the HM540 and EPIC manifests into working memory
#'
#' @param manifest_char Which manifest should be loaded? Currently, this package
#'   has been tested to work with 450k and EPIC arrays for HG19 and HG38.
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @export
#'
#' @details This function assumes that the .onLoad() function has executed 
#'   properly and (therefore) that the necessary data sets are in the cache.
#' 
#' 
#' @examples 
#'   hm450k_gr <- ImportSesameData("HM450.hg19.manifest")
#'   head(hm450k_gr)
#'

ImportSesameData <- function(manifest_char) {
  # Tiago wrote the guts of this function. See:
  # https://github.com/TransBioInfoLab/MethReg/commit/7284ce917735322762ffc3d3807d858d260c63c6
  
  # Query Data
  ehub <- ExperimentHub::ExperimentHub(localHub = FALSE)
  query <- AnnotationHub::query(ehub, c("sesameData", manifest_char))
  query <- query[query$title == manifest_char, ]
  
  # Take most recent version
  whenAdded_char <- query$rdatadateadded
  whichRecent_lgl <- whenAdded_char == max(as.Date(whenAdded_char))
  ah_id <- query$ah_id[whichRecent_lgl]
  
  # Import
  ehub[[ah_id, verbose = FALSE]]
  
}
