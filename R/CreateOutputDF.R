#' Create Output Dataframe
#'
#' @param keepCpGs_df a data frame with \code{CpG} = CpG name, \code{keep} =
#'   indicator for co-methylated CpGs, and \code{r_drop} = correlation between 
#'   the CpG with rest of the CpGs
#' @param keepContiguousCpGs_df a data frame with \code{ProbeID} = CpG name and
#'   \code{Subregion} = contiguous comethylated subregion number
#' @param CpGsOrdered_df a data frame of CpG location with \code{chr} =
#'   chromosome number, \code{pos} = genomic position, and \code{cpg} = CpG name
#' @param returnAllCpGs indicates if outputting all the CpGs in the region when
#'   there is not a contiguous comethylated region or only the CpGs in the
#'   contiguous comethylated regions
#'
#' @return a data frame with \code{CpG} = CpG name, \code{Chr} = chromosome
#'   number, \code{MAPINFO} = genomic position, \code{r_drop} = correlation
#'   between the CpG with rest of the CpGs, \code{keep} = indicator for
#'   co-methylated CpG, and \code{keep_contiguous} = contiguous comethylated
#'   subregion number
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#'    data(betasChr22_df)
#'    CpGsChr22_char <- c(
#'      "cg02953382", "cg12419862", "cg24565820", "cg04234412", "cg04824771",
#'      "cg09033563", "cg10150615", "cg18538332", "cg20007245", "cg23131131",
#'      "cg25703541"
#'    )
#'       
#'    CpGsOrdered_df <- OrderCpGsByLocation(
#'       CpGsChr22_char, arrayType="450k", output = "dataframe"
#'    )
#'    betaCluster_mat <- t(betasChr22_df[CpGsOrdered_df$cpg, ])
#'    keepCpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaCluster_mat)
#'    keepContiguousCpGs_df <- FindComethylatedRegions(CpGs_df = keepCpGs_df)
#'    CreateOutputDF(keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df)
#'    
CreateOutputDF <- function(
  keepCpGs_df,
  keepContiguousCpGs_df,
  CpGsOrdered_df,
  returnAllCpGs = FALSE
){

  if (returnAllCpGs == FALSE & all(keepContiguousCpGs_df$Subregion == 0)) {
    return(NULL)
  } 
  
  output_df <- merge(
    keepCpGs_df,
    keepContiguousCpGs_df,
    by.x = "CpG", by.y = "ProbeID",
    all.x = TRUE
  )
  output2_df <- merge(
    CpGsOrdered_df,
    output_df,
    by.x = "cpg", by.y = "CpG"
  )
  
  output3_df <- output2_df[order(output2_df$chr, output2_df$pos), ]
  output3_df[is.na(output3_df)] <- 0
  
  coMethCpGs_df <- data.frame(
    Region = NameRegion(CpGsOrdered_df),
    CpG = output3_df$cpg,
    Chr = output3_df$chr,
    MAPINFO = output3_df$pos,
    r_drop = output3_df$r_drop,
    keep = output3_df$keep,
    keep_contiguous = output3_df$Subregion
  )
  
  coMethCpGs_df

}
