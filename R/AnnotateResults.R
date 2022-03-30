#' Annotate \code{coMethDMR} Pipeline Results
#'
#' @description Given a data frame with regions in the genome, add gene symbols,
#'   UCSC reference gene accession, UCSC reference gene group and
#'   relation to CpG island.
#'
#' @param lmmRes_df A data frame returned by \code{\link{lmmTestAllRegions}}.
#' This data frame must contain the following columns:
#'    \itemize{
#'      \item{\code{chrom} : }{the chromosome the region is on, e.g. ``chr22''}
#'      \item{\code{start} : }{the region start point}
#'      \item{\code{end} : }{the region end point}
#'    }
#' @param arrayType Type of array: 450k or EPIC
#' @param nCores_int Number of computing cores to be used when executing code
#'    in parallel. Defaults to 1 (serial computing).
#' @param ... Dots for additional arguments passed to the cluster constructor.
#'    See \code{\link{CreateParallelWorkers}} for more information.
#'
#' @return A data frame with
#'    \itemize{
#'      \item the location of the genomic region's chromosome (\code{chrom}),
#'        start (\code{start}), and end (\code{end});
#'    \item UCSC annotation information (\code{UCSC_RefGene_Group},
#'        \code{UCSC_RefGene_Accession}, and \code{UCSC_RefGene_Name}); and
#'      \item a list of all of the probes in that region (\code{probes}).
#'    }
#'
#' @details The region types include \code{"NSHORE"}, \code{"NSHELF"},
#'    \code{"SSHORE"}, \code{"SSHELF"}, \code{"TSS1500"}, \code{"TSS200"},
#'    \code{"UTR5"}, \code{"EXON1"}, \code{"GENEBODY"}, \code{"UTR3"}, and
#'    \code{"ISLAND"}.
#'
#' @export
#'
#' @examples
#'    lmmResults_df <- data.frame(
#'      chrom = c("chr22", "chr22", "chr22", "chr22", "chr22"),
#'      start = c("39377790", "50987294", "19746156", "42470063", "43817258"),
#'      end   = c("39377930", "50987527", "19746368", "42470223", "43817384"),
#'      regionType = c("TSS1500", "EXON1", "ISLAND", "TSS200", "ISLAND"),
#'      stringsAsFactors = FALSE
#'    )
#'
#'    AnnotateResults(
#'      lmmRes_df = lmmResults_df,
#'      arrayType = "450k"
#'    )
#'
AnnotateResults <- function(
    lmmRes_df,
    arrayType = c("450k", "EPIC"),
    nCores_int = 1L,
    ...
){
  
    ###  Check Inputs  ###
    stopifnot(
        "data.frame" %in% class(lmmRes_df),
        all(c("chrom", "start", "end") %in% colnames(lmmRes_df))
    )
    arrayType <- match.arg(arrayType)

    lmmRes_df$start <- as.integer(lmmRes_df$start)
    lmmRes_df$end   <- as.integer(lmmRes_df$end)

    
    ###  Pull Database  ###
    switch(arrayType,
      "450k" = {
        locations_df <- 
          IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
        UCSCinfo_df <- 
          IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
        IslandsUCSCinfo_df <-
          IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC
      },
      "EPIC" = {
        locations_df <- 
          IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
        UCSCinfo_df <- 
          IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other
        IslandsUCSCinfo_df <- 
          IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC
      }
    )

    # Locations
    locations_df <- as.data.frame(locations_df)
    locations_df$cpg <- row.names(locations_df)
    rownames(locations_df) <- NULL

    # UCSC Gene Info
    UCSCinfo_df <- as.data.frame(UCSCinfo_df)
    interestingColumns_char <- c(
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Accession",
        "UCSC_RefGene_Group"
    )
    UCSCinfo_df <- UCSCinfo_df[, interestingColumns_char]

    # UCSC Island Info
    IslandsUCSCinfo_df <- as.data.frame(IslandsUCSCinfo_df)

    
    ###  Work and Return  ###
    cluster <- CreateParallelWorkers(nCores_int, ...)

    resultsAnno_ls <- bplapply(
      seq_len(nrow(lmmRes_df)),
      function(row){
        .AnnotateRow(
          row_df = lmmRes_df[row, ],
          loc_df = locations_df,
          info_df = UCSCinfo_df,
          island_df = IslandsUCSCinfo_df
        )
      },
      BPPARAM = cluster
    )
    
    do.call(rbind, resultsAnno_ls)
    
}

.AnnotateRow <- function(row_df, loc_df, info_df, island_df){
  # browser()
  
  
  ###  Filter Data Frames  ###
  # Extract Row Region
  chr   <- row_df$chrom
  start <- row_df$start
  end   <- row_df$end
  
  # Find Probes in that Region
  chr_df  <- loc_df[loc_df$chr == chr, ]
  inRegion_lgl <- chr_df$pos >= start & chr_df$pos <= end
  out_df <- chr_df[inRegion_lgl, ]
  probes_char <- out_df$cpg
  
  # Find UCSC Annotation Information for those Probes
  infoOut_df <- info_df[probes_char, ]
  
  # Find UCSC Relation to Island Information for those Probes
  islandOut_df <- island_df[probes_char, ]
  
  
  ###  Wrangle UCSC Annotation  ###
  refGeneGroup_char <- .ExtractUCSCinfo(infoOut_df$UCSC_RefGene_Group)
  refGeneAcc_char   <- .ExtractUCSCinfo(infoOut_df$UCSC_RefGene_Accession)
  refGeneName_char  <- .ExtractUCSCinfo(infoOut_df$UCSC_RefGene_Name)
  
  refIslandRelation_char <- sort(unique(islandOut_df$Relation_to_Island))
  
  
  ###  Return Annotated 1-Row Data Frame  ###
  row_df$UCSC_RefGene_Group <-
    paste0(unique(refGeneGroup_char), collapse = ";")
  row_df$UCSC_RefGene_Accession <-
    paste0(unique(refGeneAcc_char), collapse = ";")
  row_df$UCSC_RefGene_Name <-
    paste0(unique(refGeneName_char), collapse = ";")
  row_df$Relation_to_Island <-
    paste0(unique(refIslandRelation_char), collapse = ";")
  
  row_df
  
}

.ExtractUCSCinfo <- function(infoCol) {
  sort(
    unique(
      unlist(
        strsplit(infoCol, ";")
      )
    )
  )
}
