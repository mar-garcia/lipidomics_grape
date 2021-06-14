# Target ----
target <- readxl::read_xlsx("../target/data/data_tissues.xlsx", sheet = "raw_data")
colnames(target)[1] <- "ID"
target <- target[,c(grep("ID|entire", colnames(target)))]

## Compound data ----
#cmps <- read.csv("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/compounds_target.csv")
cmps <- read.csv("compounds_target.csv")

#fml.list <- list()
#for(i in seq(nrow(cmps))){
#  fml.list[[i]] <- countElements(cmps$formula[i])
#}
#tmp <- bind_rows(fml.list)
#tmp[is.na(tmp)] <- 0


# Compound ----
i <- 35
target$ID[order(-rowMeans(target[,-1]))[i]]
cmps[cmps$ID == target$ID[order(-rowMeans(target[,-1]))[i]],]

feat$rtmed <- feat$rtmed/60
feat[unlist(matchWithPpm(802.5604, feat$mzmed, ppm = 5)),]
feat_grape$rtmed[unlist(matchWithPpm(802.5604, feat_grape$mzmed, ppm = 5))]/60


# Untarget ----
feat <- feat[order(-feat$mean),]
feat[which(is.na(feat$compound))[1], ]

cmps[cmps$class == "PC" & cmps$formula == "C42H80NO8P",]
