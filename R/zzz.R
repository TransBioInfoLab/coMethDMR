#' @title Cache sesameData at Package Load
#' 
#' @description Check if the user has both the HM540 and EPIC manifests in their
#'   cache. The contents of the cache are checked via a call to the
#'   \code{ExperimentHub} function. If not all data components are available in
#'   the cache for these two platforms, we call \code{sesameDataCache} to save
#'   them to the cache for later use.
#'
#' @param libname path to package library
#' @param pkgname package name
#' 
#' @importFrom sesameData sesameDataCache sesameDataList
#' @importFrom ExperimentHub ExperimentHub
#'
#' @details arguments are unused
#' 
#' @keywords internal
#' @name coMethDMR_setup
#'
NULL

.onLoad <- function(libname, pkgname){
  
  ###  Check if data has been cached  ###
  packageStartupMessage("Checking for cached SeSAMe data.")
  platformInfo_df <- sesameData::sesameDataList()
  cachedDataNames_char <- names(ExperimentHub::ExperimentHub(localHub = TRUE))
  
  # 450k
  hmTitles_char <- c("HM450.hg19.manifest", "HM450.hg38.manifest")
  ehNames450k_char <- platformInfo_df[
    platformInfo_df$"Title" %in% hmTitles_char, "EHID", drop = TRUE
  ]
  all450k_lgl <- all(ehNames450k_char %in% cachedDataNames_char)
  
  # EPIC
  epicTitles_char <- c("EPIC.hg19.manifest", "EPIC.hg38.manifest")
  ehNamesEPIC_char <- platformInfo_df[
    platformInfo_df$"Title" %in% epicTitles_char, "EHID", drop = TRUE
  ]
  allEPIC_lgl <- all(ehNamesEPIC_char %in% cachedDataNames_char)
  
  
  ###  Download manifests as necessary  ###
  if(!all450k_lgl) {
    packageStartupMessage("Caching SeSAMe data for 450k arrays.")
    sesameDataCache("HM450")
  }
  if(!allEPIC_lgl) {
    packageStartupMessage("Caching SeSAMe data for EPIC arrays.")
    sesameDataCache("EPIC")
  }
  
  invisible()
}