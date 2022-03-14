setwd("~/GitHub/lipidomics_grape")
library(Spectra)
load("data/MS2comb.RData")
label_fun <- function(x) {
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 7)
  mzs[ints < 0.1] <- ""
  mzs
}

# normalize intensities
norm_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}
sps <- addProcessing(spsx, norm_int)

# filter low-intensity peaks
low_int <- function(x, ...) {
  x > max(x, na.rm = TRUE) * 0.01
}
sps <- filterIntensity(sps, intensity = low_int)

# get the mz of the maximum peak from each spectrum:
max_peak <- function(x, ...) {
  x[which.max(x[, 2]), , drop = FALSE]
}
sps_2 <- addProcessing(sps, max_peak)
#(idx <- which(!MsCoreUtils::between(unlist(mz(sps_2)), 397.3835  + 0.01 * c(-1, 1))))

# calculate neutral loss spectra
neutral_loss <- function(x, spectrumMsLevel, precursorMz, ...) {
  if (spectrumMsLevel == 2L) {
    x[, "mz"] <- precursorMz - x[, "mz"]
    x <- x[order(x[, "mz"]), , drop = FALSE]
  }
  x
}
sps_nl <- addProcessing(sps, neutral_loss,
                       spectraVariables = c("msLevel", "precursorMz"))

# Select spectra of interest
myidx <- which(sps$compound == "mPA 38:4" & sps$polarity == 1L)
mysps <- sps[myidx]
plotSpectra(mysps, labels = label_fun, labelPos = 2,
            labelOffset = 0.2, labelSrt = -30, labelCex = 0.6)

# Compare spectra X with the other ones
cormat <- compareSpectra(mysps, sps, ppm = 40)
idx <- order(cormat, decreasing = TRUE)[3]
sps$compound[idx]
plotSpectraMirror(mysps, sps[idx],
                  ppm = 20, labels = label_fun, labelPos = 2,
                  labelOffset = 0.2, labelSrt = -30, labelCex = 0.7)
grid()
tmp <- data.frame(c = cormat, cmp = spsx$compound)
tmp$c <- as.numeric(tmp$c)
head(tmp[order(tmp$c, decreasing = TRUE),])

cormat_nl <- compareSpectra(sps_nl[myidx], sps_nl, ppm = 40)
idx_nl <- order(cormat_nl, decreasing = TRUE)[2]
spsx$compound[idx_nl]
plotSpectraMirror(sps_nl[myidx], sps_nl[idx_nl],
                  ppm = 20, labels = label_fun, labelPos = 2,
                  labelOffset = 0.2, labelSrt = -30, labelCex = 0.7)
grid()
tmp <- data.frame(c = cormat_nl, cmp = spsx$compound)
tmp$c <- as.numeric(tmp$c)
head(tmp[order(tmp$c, decreasing = TRUE),])


library(MsCoreUtils)
ions[MsCoreUtils::between(ions$mz, 255.2327 + 0.01 * c(-1,1)) & nl$positive == FALSE,]
nl$massx <- (nl$mass_add + 1.007276)*(-1)
nl[MsCoreUtils::between(nl$massx, 414.2181 + 0.01 * c(-1,1)) & nl$positive == FALSE,]
