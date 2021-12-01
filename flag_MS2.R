library(tidyverse)
library(Spectra)
plotms2 <- function(mz, int){
  plot(mz, int, type = "h")
  idx <- which(int / max(int) > 0.1)
  text(mz[idx], int[idx], round(mz, 4)[idx], cex=0.8)
}
.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

s <- "tissues"
p <- "POS"
cmp <- "PE_16:0_18:2 _I"
prec <- 716.5180
rt <- 16.62
frag <- 575.503387


load(paste0("data/RData/data_", s, "_", p, ".RData"))
ms2_proc <- ms2
load(paste0("data/RData/MS2_", s, "_", p, ".RData"))
ms2 <- ms2[!grepl("STDmix", basename(ms2@backend@spectraData$dataOrigin))]

i <- which(db_exp$compound == cmp)[1] 
ad <- db_exp$adduct[i]
ms2sub <- ms2_proc[grep(db_exp$FT[i], ms2_proc$FT)]
#mypairs <- t(combn(seq(length(ms2sub)), 2))
#pippo <- mypairs %>%
#  as_tibble() %>% 
#  mutate(mylist = map2(V1,V2, function(x,y) c(x,y))) %>% pull(mylist)
#pluto <- lapply(pippo, function(x) {
#  list(ms2sub[x[1]], ms2sub[x[2]])
#}) 
#myf <- function(x){
#  Spectra::compareSpectra(x[[1]], x[[2]], ppm = 10)
#}
#res <- unlist(lapply(pluto, myf))
#hist(res)
ms2comb <- Spectra::combineSpectra(
  ms2sub, intensityFun = base::sum, mzFun = base::mean, 
  tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
  weighted = TRUE)
dt <- as.data.frame(cbind("mz" = unlist(mz(ms2comb)),
                          "intensity" = unlist(intensity(ms2comb))))
plotms2(dt$mz, dt$intensity)
j <- 1
dt <- as.data.frame(cbind("mz" = unlist(mz(ms2sub[j])),
                          "intensity" = unlist(intensity(ms2sub[j]))))
plotms2(dt$mz, dt$intensity)

ms2sub <- filterPrecursorMz(ms2, prec + 0.01 * c(-1, 1))
ms2sub <- filterRt(ms2sub, rt*60 + 10 * c(-1, 1))
table(round(ms2sub$mz_max))
#idx <- which(ms2sub$mz_max > (floor(frag)+1) | ms2sub$mz_max < floor(frag))
idx <- unlist(matchWithPpm(frag, ms2sub$mz_max, ppm = 10))
idx <- seq(length(ms2sub))[!seq(length(ms2sub)) %in% idx]
# TAGs ----
idx <- unlist(matchWithPpm(c(601.5195, 575.5034), ms2sub$mz_max, ppm = 10))
table(round(ms2sub$mz_max[idx]))
idx <- seq(length(ms2sub))[!seq(length(ms2sub)) %in% idx]
#length(ms2sub)
# -----
length(idx)
write.csv(data.frame(
  "File" = basename(ms2sub@backend@spectraData@listData$dataOrigin)[idx],
  "scan" = ms2sub@backend@spectraData@listData$scanIndex[idx],
  "compound" = rep(cmp, length(idx)),
  "adduct" = rep(ad, length(idx)),
  "observations" = rep("X", length(idx)),
    #paste0("Main ion should be ", sprintf("%.5f", frag), 
     #      " instead of ", sprintf("%.5f", ms2sub$mz_max[idx])),
  "counts" = ms2sub$intensity_max[idx],
  "polarity" = rep(p, length(idx)),
  "study" = rep(s, length(idx))
), "x.csv")


ms2subx <- ms2sub
table(round(ms2subx$mz_max))
ms2sub <- ms2subx[round(ms2subx$mz_max) == 602]
ms2comb <- Spectra::combineSpectra(
  ms2sub, intensityFun = base::sum, mzFun = base::mean, 
  tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
  weighted = TRUE)
dt <- as.data.frame(cbind("mz" = unlist(mz(ms2comb)),
                          "intensity" = unlist(intensity(ms2comb))))
plotms2(dt$mz, dt$intensity)

