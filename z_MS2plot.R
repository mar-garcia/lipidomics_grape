library(MsCoreUtils)
plotms2 <- function(mz, int, main = main){
  plot(mz, int, main = main, type = "h")
  idx <- which(int / max(int) > 0.1)
  text(mz[idx], int[idx], round(mz, 4)[idx], cex=0.8)
}
max_peak <- function(x, ...) {
  x[which.max(x[, 2]), , drop = FALSE]
}
load("data/RData.RData")


s <- "maturation" # specify "maturation" or "tissues"
p <- "POS" # specify "POS" or "NEG"
dda_spectra <- get(paste("dda_spectra", s, p, sep = "_"))
dda_xdata <- get(paste("dda_xdata", s, p, sep = "_"))

ft_mz <- 782.6438    
ft <- get(paste("ft", s, p, sep = "_"))
ft_id <- rownames(featureDefinitions(dda_xdata, mz = ft_mz, ppm = 10))
cp_id <- rownames(chromPeaks(dda_xdata)[ft[ft_id, "peakidx"][[1]],])
(ft_spectra <- dda_spectra[dda_spectra$peak_id %in% cp_id])
dev.off()
for(j in rev(seq(length(ft_spectra)))){
  dt <- as.data.frame(cbind("mz" = unlist(mz(ft_spectra[j])),
                            "intensity" = unlist(intensity(ft_spectra[j]))))
  plotms2(dt$mz, dt$intensity, main = j)
}

#### STDS #####
load("data/STDmix_MS2.RData")
(i_spectra <- filterPrecursorMz(ms2, ft_mz + 0.01 * c(-1, 1)))
(ft_spectra <- filterRt(i_spectra,21.81*60 + 10 * c(-1, 1)))


# define a function that returns the maximum peak from each spectrum:
sps_2 <- addProcessing(ft_spectra, max_peak)
idx <- which(!between(unlist(mz(sps_2)), 603.534687 + 0.01 * c(-1, 1)))


write.csv(data.frame(
  "File" = basename(ft_spectra@backend@spectraData@listData$dataOrigin)[idx],
  "scan" = ft_spectra@backend@spectraData@listData$scanIndex[idx],
  "observations" =  "precursor is a (2)13C ion"
  #"noisy MS2"
  #"intensity of 13C[M+H]+ very high"
  #"main ion (760.5753) is higher than precursor (758.56943)"
    #paste("Main peak should be", 603.5347, "instead of", round(unlist(mz(sps_2))[idx], 4))
), "x.csv")
