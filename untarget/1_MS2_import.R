library(Spectra)
s <- "maturation"
#for(s in c("tissues", "maturation")){
p <- "POS"
#for(p in c("NEG", "POS")){
fls <- list.files(paste0(s, "/data/", p, "_DDA_mzML/"), full.names = TRUE)
ms2 <- Spectra(fls, backend = MsBackendDataFrame())
ms2 <- ms2[msLevel(ms2) == 2]
mz_max <- rep(NA, length(ms2))
intensity_max <- rep(NA, length(ms2))
startpoint <- Sys.time()
for(i in 1:length(ms2)){
  dt <- cbind(unlist(mz(ms2[i])), unlist(intensity(ms2[i])))
  idx <- which.max(dt[,2])
  mz_max[i] <- dt[idx, 1]
  intensity_max[i] <- dt[idx, 2]
  print(i)
}
Sys.time()-startpoint
ms2$mz_max <- mz_max
ms2$intensity_max <- intensity_max
save(ms2, file = paste0(s, "/data/RData/MS2_", p, ".RData"))




