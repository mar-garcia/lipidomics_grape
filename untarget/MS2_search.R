library(Spectra)
library(xcms)
library(MsCoreUtils)
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


ms2sub <- filterPrecursorMz(ms2_POS, 871.5715 + 0.01 * c(-1, 1))
ms2sub <- filterRt(ms2sub, 18.98*60 + 10 * c(-1, 1))
#ms2sub <- ms2sub[containsMz(ms2sub,c(570.5049), tolerance = 0.005)]
#ms2sub <- ms2sub[containsMz(ms2sub, c(491.3220), tolerance = 0.005)]
length(ms2sub)
intensitats <- c()
for(i in seq(length(ms2sub))){
  idx <- c(
    which(unlist(mz(ms2sub[i])) > 872),
    which(unlist(mz(ms2sub[i])) > 283.5 & unlist(mz(ms2sub[i])) < 283.9),
    which(unlist(mz(ms2sub[i])) > 341.1 & unlist(mz(ms2sub[i])) < 341.3)
  )
  int.noise <- max(unlist(intensity(ms2sub[i]))[idx])
  int.good <- max(unlist(intensity(ms2sub[i]))[-idx])
  intensitats <- c(intensitats, int.good / int.noise)
}
rm(tmp)
tmp <- order(intensitats)[(length(intensitats)-10):length(intensitats)]
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
  
  xdata <- readMSData(ms2sub[j]@backend@spectraData$dataOrigin, mode = "onDisk")
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
                    sprintf("%.2f", ms2sub[j]@backend@spectraData$rtime/60)),
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
tmp2[order(tmp2)]
