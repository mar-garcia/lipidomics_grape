library(MsCoreUtils)
library(Spectra)
.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}
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
ft <- get(paste("ft", s, p, sep = "_"))

#ft_mz <- 981.5788146
#ft_rt <- 14.74978333
#ft_id <- "FT0580"
#(ft_id <- rownames(featureDefinitions(dda_xdata, mz = ft_mz, ppm = 10)))
#(ft_id <- pks$FT[pks$dataset == paste(s, p, sep = "_") & MsCoreUtils::between(pks$mzmed, ft_mz + 0.01 * c(-1,1)) &
#                   MsCoreUtils::between(pks$rtmed, ft_rt + 10/60 * c(-1, 1))])
cp_id <- rownames(chromPeaks(dda_xdata)[ft[ft_id, "peakidx"][[1]],])
(ft_spectra <- dda_spectra[dda_spectra$peak_id %in% cp_id])
dev.off()
for(j in rev(seq(length(ft_spectra)))){
  dt <- as.data.frame(cbind("mz" = unlist(mz(ft_spectra[j])),
                            "intensity" = unlist(intensity(ft_spectra[j]))))
  plotms2(dt$mz, dt$intensity, main = j)
  #abline(h = max(unlist(intensity(ft_spectra[j])))/2)
}

#### STDS #####
load("data/STDmix_MS2.RData")
(i_spectra <- filterPrecursorMz(ms2, ft_mz + 0.01 * c(-1, 1)))
(ft_spectra <- filterRt(i_spectra,ft_rt*60 + 10 * c(-1, 1)))


# define a function that returns the maximum peak from each spectrum:
sps_2 <- addProcessing(ft_spectra, max_peak)
(idx <- which(!MsCoreUtils::between(unlist(mz(sps_2)), 506.325213+1.006277*7 + 0.01 * c(-1, 1))))

l <- matchWithPpm(unlist(mz(sps_2)), c(267.232954, 403.225499), ppm = 10)
idx <- which(lapply(l,length) == 0)

write.csv(data.frame(
  "File" = basename(ft_spectra@backend@spectraData@listData$dataOrigin)[idx],
  "scan" = ft_spectra@backend@spectraData@listData$scanIndex[idx],
  "observations" =  
    "precursor is a (2)13C ion"
    #"noisy MS2"
    #"outlier MS2"
    #"intensity of 13C very high"
  #"main ion (781.5597) is higher than precursor (780.5524)"
  #paste("Main peak should be", 513.3692, "instead of", round(unlist(mz(sps_2))[idx], 4))
), "x.csv")




ms2comb <- combineSpectra(
  ft_spectra[idx], intensityFun = base::sum, mzFun = base::mean, 
  tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
  weighted = TRUE)
if(length(ms2comb) > 1){
  ms2comb <- combineSpectra(
    ms2comb, intensityFun = base::sum, mzFun = base::mean, 
    tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
    weighted = TRUE)
}
dt <- as.data.frame(cbind("mz" = unlist(mz(ms2comb)),
                          "intensity" = unlist(intensity(ms2comb))))
plotms2(dt$mz, dt$intensity, main ="X")




norm_intensities <- function(x, ...) {
  x[, 2] <- x[, 2] / max(x[,2])
  x
}
sps_2 <- addProcessing(ft_spectra , norm_intensities)
sps_2 <- replaceIntensitiesBelow(sps_2, threshold = 0.1, value = 0)
sps_2 <- filterIntensity(sps_2, intensity = c(0.1, Inf))
l <- matchWithPpm(mz(sps_2), c(502.329201, 520.339765), ppm = 10)
(idx <- which(lapply(l,length) == 0))
