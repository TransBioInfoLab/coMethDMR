
islands <- readRDS("ISLAND_3_200.rds")

islandName <- names(islands)
islandNameChr22 <- islandName[grep("chr22:", islandName)][1:20]
islandsChr22 <- islands[islandNameChr22]
saveRDS(islandsChr22, "CpGislandsChr22_ex.RDS")
