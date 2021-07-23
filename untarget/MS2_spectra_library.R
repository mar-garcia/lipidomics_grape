library(Spectra)
norm_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}

load("maturation/data/RData/MS2_library_POS.RData")
ms2_POS <- ms2
load("tissues/data/RData/MS2_library_POS.RData")
ms2_POS <- c(ms2_POS, ms2)
load("maturation/data/RData/MS2_library_NEG.RData")
ms2_NEG <- ms2
load("tissues/data/RData/MS2_library_NEG.RData")
ms2_NEG <- c(ms2_NEG, ms2)
rm(ms2)

ms2_dt <- openxlsx::read.xlsx("output/MS2_spectra.xlsx")
ms2_dt <- ms2_dt[!grepl("Unk|TAG|DAG", ms2_dt$compound),]
cmps <- unique(ms2_dt$compound)

for(j in 1:2){
  p <- c("POS", "NEG")[j]
  i.ms2 <- get(paste0("ms2_", p))
  for(i in seq(length(cmps))){
    i.ms2.dt <- ms2_dt[ms2_dt$compound == cmps[i] & ms2_dt$polarity == p,]
    if(nrow(i.ms2.dt) > 0){
      i.add <- unique(i.ms2.dt$adduct)
      for(k in seq(length(i.add))){
        i.ms2.dt.k <- i.ms2.dt[i.ms2.dt$adduct == i.add[k], ]
        i.ms2x <- i.ms2[paste(basename(dataOrigin(i.ms2)), scanIndex(i.ms2)) %in% 
                          paste(i.ms2.dt.k$file, i.ms2.dt.k$scan)]
        i.ms2x <- addProcessing(i.ms2x, norm_int)
        i.ms2x <- Spectra::combineSpectra(
          i.ms2x, intensityFun = base::mean, mzFun = base::mean, 
          tolerance = 0.01, minProp = 0.5, peaks = "intersect", 
          weighted = TRUE)
        i.ms2x$name <- cmps[i]
        i.ms2x$adduct <- i.add[k]
        if(exists(paste0("sps_ms2_", p))){
          sps_ms2 <- get(paste0("sps_ms2_", p))
          sps_ms2 <- c(sps_ms2, i.ms2x)
          assign(paste0("sps_ms2_", p), sps_ms2)
        } else {
          sps_ms2 <- i.ms2x
          assign(paste0("sps_ms2_", p), sps_ms2)
        }
      }
    }
  }
}
save(sps_ms2_NEG, sps_ms2_POS, file = "output/MS2_spectra.RData")
save(sps_ms2_NEG, sps_ms2_POS, file = "../../lipidomics_tool/MS2_spectra.RData")
