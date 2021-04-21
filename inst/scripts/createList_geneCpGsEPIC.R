library(coMethDMR)
library(tidyr)

annotDF <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other

cpgGene_df <- subset(annotDF, select="UCSC_RefGene_Name")

cpgGene_df <- as.data.frame(subset(cpgGene_df, UCSC_RefGene_Name != ""))

cpgGeneCG_df <- subset(cpgGene_df, substr(row.names(cpgGene_df),1,2) == "cg")

cpgGeneCG_df$CpG <- rownames(cpgGeneCG_df)

cpgGeneCGlong_df <- unique(separate_rows(cpgGeneCG_df,1,sep=";"))

AllGeneNames <- unique(cpgGeneCGlong_df$UCSC_RefGene_Name)

### make a list, where each item include cpgs for each gene

allGeneRegions_ls <- list()
for (i in 1:length(AllGeneNames)){

  cpg_df <- subset(cpgGeneCGlong_df, UCSC_RefGene_Name == AllGeneNames[i])

  geneRegion <- cpg_df$CpG

  geneRegion_ls <- list(geneRegion)
  names(geneRegion_ls) <- AllGeneNames[i]

  allGeneRegions_ls <- c(allGeneRegions_ls, geneRegion_ls)

}

region3 <- allGeneRegions_ls [lapply(allGeneRegions_ls, length) >=3]
saveRDS(region3, "inst/extdata/EPIC_GeneByName_3.rds")

region3_200 <- lapply(region3,
                      CloseBySingleRegion,
                      arrayType = "EPIC",
                      maxGap = 200,
                      minCpGs = 3)

region3_200 <- unlist(region3_200, recursive=FALSE)

region3_200Unique <- unique(region3_200)

region3_200_df <- lapply(region3_200Unique,
                         OrderCpGsByLocation,
                         arrayType = c("EPIC"),
                         output = "dataframe")

names(region3_200Unique) <- lapply(region3_200_df, NameRegion)

saveRDS(region3_200Unique, "Gene_3_200EPIC.rds")

### remove chrx and chrY ###

region <- readRDS("Gene_3_200EPIC.rds")

regionNames <- names(region)

chrX_Y <- grep(paste("chrX", "chrY", sep="|"), regionNames)

regionNoXY <- region[-chrX_Y]

saveRDS(regionNoXY, "Gene_3_200EPICnoXY.rds")


