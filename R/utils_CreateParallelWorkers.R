#' Create a Parallel Computing Cluster
#'
#' @description This function is an operating-system agnostic wrapper for the
#'   \code{\link[BiocParallel]{SnowParam}} and
#'    \code{\link[BiocParallel]{MulticoreParam}} constructor functions.
#'
#' @param nCores The number of computing cores to make available for coMethDMR
#'   computation
#' @param ... Additional arguments passed to the cluster constructors.
#' 
#' @return A parameter class for use in parallel evaluation
#'
#' @details This function checks the operating system and then creates a
#'   cluster of workers using the \code{SnowParam} function for Windows machines
#'   and the \code{MulticoreParam} function for non-Windows machines.
#'
#' @export
#'
#' @importFrom BiocParallel SnowParam SerialParam MulticoreParam
#' @examples
#'    workers_cl <- CreateParallelWorkers(nCores = 4)
#'    
CreateParallelWorkers <- function(nCores, ...){

  if (nCores > 1L) {
    if (Sys.info()["sysname"] == "Windows"){
      SnowParam(workers = nCores, type = "SOCK", ...)
    } else {
      MulticoreParam(workers = nCores, ...)
    }
  } else {
    SerialParam(...)
  }

}
