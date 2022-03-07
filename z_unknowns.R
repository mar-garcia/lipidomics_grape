library(CompoundDb)
library(MetaboAnnotation)
cmps <- compound_tbl_sdf("data/LMSD_structures.sdf")
write.csv(as.data.frame(cmps[,c(1:6,8)]), "data/LMSD_structures.csv")


cmps <- read.table("data/LMSD_COMP_DB_DATA.tsv", sep = "\t", header = TRUE)
idx <- grep("mass", colnames(cmps))
colnames(cmps)[idx] <- "exactmass"

cmps <- read.csv("data/LMSD_structures.csv")


unks <- c(658.5387)
param <- Mass2MzParam(ppm = 10, adducts = c("[M+H]+", "[M+NH4]+", "[M+K]+"))
(pks <- as.data.frame(matchedData(matchMz(unks, cmps, param = param))))[,c(1:4,10)]

load("data/ionsdb.RData")
param2 <- MzParam(ppm = 10)
(pks2 <- as.data.frame(matchedData(matchMz(unks, ionsdb, param = param2,
                                           mzColname = c("ion_mz")))))

