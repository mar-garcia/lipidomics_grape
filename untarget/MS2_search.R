library(Spectra)
library(xcms)
library(MsCoreUtils)
.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}
inj <- readxl::read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/standards_injections.xlsx")
inj <- inj[!is.na(inj$ID_cmp),]

load("maturation/data/RData/MS2_library_POS.RData")
ms2_POS <- ms2
load("tissues/data/RData/MS2_library_POS.RData")
ms2_POS <- c(ms2_POS, ms2)
load("maturation/data/RData/MS2_library_NEG.RData")
ms2_NEG <- ms2
load("tissues/data/RData/MS2_library_NEG.RData")
ms2_NEG <- c(ms2_NEG, ms2)
rm(ms2)


ms2sub <- filterPrecursorMz(ms2_NEG, 447.348 + 0.005 * c(-1, 1))
ms2sub <- filterRt(ms2sub, 15.39*60 + 10 * c(-1, 1))
#ms2sub <- filterMzValues(ms2sub, 283.7, ppm = 20, keep = FALSE)

tmp <- intensity(ms2sub)[matchWithPpm(311.2956, mz(ms2sub), ppm = 10)[[1]]]/max(intensity(ms2sub))
tmp <- lapply(tmp, function(x) if (length(x) == 0) {0} else {x})
ms2sub <- ms2sub[unlist(tmp==1)]

#tmp <- mz(ms2sub)[order(intensity(ms2sub), decreasing = T)]
#tmp <- sapply(tmp, "[[", 3)
#ms2sub <- ms2sub[unlist(matchWithPpm(255.23295, tmp, ppm = 10))]

length(ms2sub)
intensitats <- c()
for(i in seq(length(ms2sub))){
  idx <- c(
    #which(unlist(mz(ms2sub[i])) > 600.4 & unlist(mz(ms2sub[i])) < 600.6),
    #which(unlist(mz(ms2sub[i])) > 575.4 & unlist(mz(ms2sub[i])) < 575.6),
    which(unlist(mz(ms2sub[i])) > mean(precursorMz(ms2sub))+1),
    which(unlist(mz(ms2sub[i])) > 283.5 & unlist(mz(ms2sub[i])) < 283.9),
    which(unlist(mz(ms2sub[i])) > 341.1 & unlist(mz(ms2sub[i])) < 341.8)
  )
  int.noise <- max(unlist(intensity(ms2sub[i]))[idx])
  int.good <- max(unlist(intensity(ms2sub[i]))[-idx])
  intensitats <- c(intensitats, int.good / int.noise)
}
rm(tmp)
tmp <- order(intensitats)[(length(intensitats)-15):length(intensitats)]
dev.off()
if(exists("tmp")){
  i.seq <- tmp
} else {
  i.seq <- seq(length(ms2sub))
}
tmp2 <- c()
par(mfrow = c(1, 2))
for(i in i.seq){
  if(exists("tmp")){
    j <- i
  } else {
    j <- order(intensitats)[i]
    #j <- order(rtime(ms2sub))[i]
  }
  
  xdata <- readMSData(#gsub("garciaalom", "lenovo", 
    ms2sub[j]@backend@spectraData$dataOrigin#)
    , mode = "onDisk")
  chr <- chromatogram(
    xdata, mz = ms2sub[j]@backend@spectraData$precursorMz + 0.01 * c(-1, 1),
    rt = ms2sub[j]@backend@spectraData$rtime + 30 * c(-1, 1)
  )
  plot(chr, xaxt="n",
       main = gsub("_DDA.mzML", "", 
                   basename(ms2sub[j]@backend@spectraData$dataOrigin)))
  axis(1, at = seq(0, 60*30, 6), labels = sprintf("%.2f", seq(0, 30, 6/60))) 
  points(ms2sub[j]@backend@spectraData$rtime, 
         intensity(chr[[1]])[closest(ms2sub[j]@backend@spectraData$rtime, 
                                     rtime(chr[[1]]))], pch = 8)
  
  plot(unlist(mz(ms2sub[j])), 
       unlist(intensity(ms2sub[j])) / max(unlist(intensity(ms2sub[j]))), 
       type = "h", 
       xlab = "mz", ylab = "rel. intensity",
       main = paste(sprintf("%.4f", ms2sub[j]@backend@spectraData$precursorMz), 
                    "@", 
                    sprintf("%.2f", ms2sub[j]@backend@spectraData$rtime/60), "-",
                    scanIndex(ms2sub[j])),
       #sub = paste0("MS", ms2sub[j]@backend@spectraData$msLevel)
  )
  idx <- which(unlist(intensity(ms2sub[j])) / max(unlist(intensity(ms2sub[j]))) > 0.1)
  text(unlist(mz(ms2sub[j]))[idx], 
       (unlist(intensity(ms2sub[j])) / max(unlist(intensity(ms2sub[j]))))[idx],
       round(unlist(mz(ms2sub[j]))[idx], 4))
  tmp2 <- c(tmp2, paste(
    inj$ID_file[
      inj$filename == gsub(".mzML", "", basename(ms2sub[j]@backend@spectraData$dataOrigin))],
    basename(ms2sub[j]@backend@spectraData$dataOrigin), 
    ms2sub[j]@backend@spectraData$acquisitionNum))
}
#tmp2[order(tmp2)]
tmp2
par(mfrow=c(1,1))

##############################################################
dt <- data.frame(cbind(unlist(mz(ms2sub[j])), 
            unlist(intensity(ms2sub[j])) / max(unlist(intensity(ms2sub[j])))))
dt$X2 <- dt$X2*100
dt <- dt[dt$X2 > 1,]
write.table(dt, "dt.txt", row.names = F, col.names = F)
##############################################################
library(pheatmap)
label_fun <- function(x) {
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 4)
  mzs[ints < 5] <- ""
  mzs
}

cormat <- compareSpectra(i.ms2, ppm = 20)
hm <- pheatmap(cormat)
for(i in hm$tree_row$order){
  plot(unlist(mz(i.ms2[i])), unlist(intensity(i.ms2[i])), type = "h", 
       #main = paste(i, round(rtime(i.ms2[i])/60, 2)), 
       main = i, xlab = "m/z", ylab = "intensity")
  text(unlist(mz(i.ms2[i])), unlist(intensity(i.ms2[i])),
       round(unlist(mz(i.ms2[i])), 4), cex = 0.8)
}
j <- 6
for(i in hm$tree_row$order){
  plotSpectraMirror(i.ms2[j], i.ms2[i], labels = label_fun, 
                    labelPos = 2, labelOffset = 0.2, labelSrt = -30)
}

i.ms2x <- combineSpectra(
  i.ms2, intensityFun = base::mean, mzFun = base::mean, 
  tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
  weighted = TRUE)
plotSpectraMirror(i.ms2[j], i.ms2x, labels = label_fun, 
                  labelPos = 2, labelOffset = 0.2, labelSrt = -30)

