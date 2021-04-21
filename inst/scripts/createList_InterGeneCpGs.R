library(coMethDMR)
library(tidyr)

### Find intergenic regions for each chr in the annotation file

annotDF <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other

cpgGene_df <- subset(annotDF, select="UCSC_RefGene_Name")

cpgInterGene_df <- as.data.frame(subset(cpgGene_df, UCSC_RefGene_Name == ""))

cpgInterGeneCG_df <- subset(cpgInterGene_df, substr(row.names(cpgInterGene_df),1,2) == "cg")

cpgInterGeneCG_df$CpG <- rownames(cpgInterGeneCG_df)

loc <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations

cpgInterGeneCGLoc_df <- merge(cpgInterGeneCG_df, loc, by.x = "CpG", by.y = "row.names")

cpgInterGeneCGLoc_df <- as.data.frame(cpgInterGeneCGLoc_df)

AllChrs <- unique(cpgInterGeneCGLoc_df$chr)

### make a list, where each item include cpgs for each chr

allChrRegions_ls <- list()
for (i in 1:length(AllChrs)){

  cpg_df <- subset(cpgInterGeneCGLoc_df, chr == AllChrs[i])

  chrRegion <- cpg_df$CpG

  chrRegion_ls <- list(chrRegion)
  names(chrRegion_ls) <- AllChrs[i]

  allChrRegions_ls <- c(allChrRegions_ls, chrRegion_ls)

}

### Create a list of close by intergenic regions

region3 <- allChrRegions_ls [lapply(allChrRegions_ls, length) >=3]
region3_200 <- lapply(region3,
                      CloseBySingleRegion,
                      arrayType = "450k",
                      maxGap = 200,
                      minCpGs = 3)

region3_200 <- unlist(region3_200, recursive=FALSE)

region3_200_df <- lapply(region3_200,
                         OrderCpGsByLocation,
                         arrayType = c("450k"),
                         output = "dataframe")

names(region3_200) <- lapply(region3_200_df, NameRegion)

saveRDS(region3_200, "InterGene_3_200.rds")

### remove chrx and chrY ###

region <- readRDS("InterGene_3_200.rds")

regionNames <- names(region)

chrX_Y <- grep(paste("chrX", "chrY", sep="|"), regionNames)

regionNoXY <- region[-chrX_Y]

saveRDS(regionNoXY, "InterGene_3_200.rds")

