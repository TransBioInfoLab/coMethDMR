# Create Chromosome 22 CpG Islands Example Data Set
# Gabriel Odom and Lissette Gomez
# UNKNOWN DATE
# UPDATED: 2021-05-04


######  Preliminary Data Setup  ###############################################
# Please download the IMA full 450k annotation data ("fullannotInd.rda") from:
# https://rforge.net/IMA/#sec-4
# 
# Once you have that data, extract and save the "ISLANDSInd" data object,
#   specifically the Site IDs (SIDs) as opposed to the Position IDs (PIDs). The 
#   code below assumes such an object exists in your current global environment.

islandRegions_ls <- ISLANDInd$SID[lengths(ISLANDInd$SID) >= 3]

library(parallel)
clust <- makeCluster(20L)

###  Find the Clusters  ###
system.time(
  islandRegions_min3_200bp_ls <- parLapply(
    cl = clust,
    X = islandRegions_ls,
    fun = coMethDMR::CloseBySingleRegion,
    arrayType = "450k",
    maxGap = 200,
    minCpGs = 3
  )
)
# Windows Parallel: 6.453333 min for 23,014 regions over 20 cores

islandClusters_min3_200bp_ls <- unlist(
  islandRegions_min3_200bp_ls,
  recursive = FALSE
)
  

###  Find the Cluster Regions  ###
system.time(
  islandClustersDF_ls <- parLapply(
    cl = clust,
    X = islandClusters_min3_200bp_ls,
    fun = coMethDMR::OrderCpGsByLocation,
    arrayType = "450k",
    output = "dataframe"
  )
)
# Windows Parallel: 5.032 min for 20,590 clusters over 20 cores

names(islandClusters_min3_200bp_ls) <- lapply(islandClustersDF_ls, NameRegion)
saveRDS(islandClusters_min3_200bp_ls, "inst/extdata/ISLAND_3_200_new.rds")
# The original has 613 fewer clusters (0.3Mb less), so we can include both.



######  Creating the Example Data Set  ########################################
islands <- readRDS("ISLAND_3_200.rds")

islandName <- names(islands)
islandNameChr22 <- islandName[grep("chr22:", islandName)][1:20]
islandsChr22 <- islands[islandNameChr22]
saveRDS(islandsChr22, "CpGislandsChr22_ex.rds")
