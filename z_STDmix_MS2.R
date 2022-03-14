# import MS2 of standards
s <- "tissues"
for(p in c("POS", "NEG")){
  fls <- list.files(paste0("data/", s, "/", p, "_DDA_mzML/"), full.names = TRUE)
  fls <- fls[grep("STDmix", fls)]
  library(Spectra)
  ms2 <- Spectra(fls, backend = MsBackendDataFrame())
  ms2 <- ms2[msLevel(ms2) == 2]
  assign(paste("ms2", p, sep = "_"), ms2)
}
ms2 <- c(ms2_POS, ms2_NEG)

ms2_exclude <- read.csv("data/MS2_exclude_STD.csv")
idx <- rep(NA, nrow(ms2_exclude))
for(i in seq(length(idx))){
  x <- which(
    basename(
      ms2@backend@spectraData@listData$dataOrigin
    ) == ms2_exclude$File[i] & 
      ms2@backend@spectraData@listData$scanIndex == ms2_exclude$scan[i]
  )
  if(length(x) > 0){
    idx[i] <- x
  }
}
idx <- idx[!is.na(idx)]
ms2 <- ms2[-idx]

save(ms2, file = "data/MS2_STDmix.RData")
