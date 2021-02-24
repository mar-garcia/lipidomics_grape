library(xcms)
library(CluMSID)
library(CompoundDb)
library(Rdisop)

polarity <- "POS" # specify "POS" or "NEG"
load(paste0("data/RData/data_XCMS_", polarity, ".RData"))
load(paste0("data/RData/MS2_library_", polarity, ".RData"))
data <- featureValues(xdata, method = "sum", value = "into")
data[is.na(data)] <- 0
features <- data.frame(featureDefinitions(xdata))

z.ft <- "FT0764"
which.max(data[z.ft,])
z.idx <- which(rownames(features) == z.ft)
# RT range
z.features <- features[(features$rtmed > (features$rtmed[z.idx] - 10)) & 
                         ( features$rtmed < (features$rtmed[z.idx] + 10)), ]
# intensity correlation
z.data <- data[rownames(data) %in% rownames(z.features), ]
z.features$cor_int <- cor(t(z.data), z.data[z.ft,])
z.features <- z.features[z.features$cor_int > 0.7, ]
# peak-shape correlation
z.xdata <- filterFile(xdata, which(xdata$order == "x094"))
z.chr <- chromatogram(z.xdata,
                      mz = features$mzmed[z.idx] + 0.01 * c(-1, 1),
                      rt = features$rtmed[z.idx] + 10 * c(-1,1),
                      aggregationFun = "max")
z.cor <- c()
for(i in seq(nrow(z.features))){
  z.chr2 <- chromatogram(z.xdata,
                         mz = z.features$mzmed[i] + 0.01 * c(-1, 1),
                         rt = z.features$rtmed[i] + 10 * c(-1,1),
                         aggregationFun = "max")
  z.cor <- c(z.cor, correlate(z.chr[[1]], z.chr2[[1]]))
}
z.cor[is.na(z.cor)] <- 0
z.features$cor_ps <- z.cor
z.features <- z.features[z.features$cor_ps > 0.7, ]
z.features$i <- seq(nrow(z.features))

z.features[,c("mzmed", "rtmed", "cor_int", "cor_ps", "i")]

tmp <- unlist(mass2mz(getMolecule("C33H58D7NO3")$exactmass, adduct = adducts()))
unlist(matchWithPpm(tmp, z.features$mzmed, ppm = 10))




z.xdata %>%
  filterRt(rt = features$rtmed[z.idx] + 10 * c(-1,1)) %>%
  filterMz(mz = features$mzmed[z.idx]+21.9819 + 0.01 * c(-1, 1)) %>%
  plot(type = "XIC")




chr <- chromatogram(xdata, 
                    mz = c(features$mzmin[z.idx]-0.01, features$mzmax[z.idx]+0.01),
                    rt = c(features$rtmin[z.idx]-10, features$rtmax[z.idx]+50))
plotChromPeakDensity(chr)




mz <- features$mzmed[z.idx]
rt <- features$rtmed[z.idx]
ms2sub <- getSpectrum(ms2list, "precursor", mz, mz.tol = 0.01)
ms2sub <- getSpectrum(ms2sub, "rt", rt, rt.tol = 10)
if(length(ms2sub) > 1){
  intensitats <- c()
  for(i in seq(ms2sub)){
    idx <- substring(gsub(".*\\.","", accessSpectrum(ms2sub[[i]])[,1]), 1, 1)>1
    int.noise <- accessSpectrum(ms2sub[[i]])[idx,2][which.max(accessSpectrum(ms2sub[[i]])[idx,2])]
    int.good <- accessSpectrum(ms2sub[[i]])[-idx,2][which.max(accessSpectrum(ms2sub[[i]])[-idx,2])]
    intensitats <- c(intensitats, int.good / int.noise)
  }
}
dev.off()
par(mfrow=c(1,2))
if(length(ms2sub) > 30){
  for(i in (length(ms2sub)-30):length(ms2sub)){
    j <- order(intensitats)[i]
    
    raw_data <- readMSData(
      files = paste0("data/", polarity, "_DDA_mzmL/", ms2sub[[j]]@annotation), 
      mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = mz + 0.01 * c(-1, 1), 
                        rt = rt + 20 * c(-1, 1)
    )
    plot(chr, xlim = rt + 20 * c(-1, 1)
    )
    abline(v=ms2sub[[j]]@rt)
    
    specplot(ms2sub[[j]],
             main = paste("id:", ms2sub[[j]]@id, " - ", 
                          ms2sub[[j]]@annotation))
    print(paste0(j, ": ", gsub(".*\\/", "", ms2sub[[j]]@annotation), " - ", ms2sub[[j]]@id, 
                 " - ", ms2sub[[j]]@rt))
  }
} else if(length(ms2sub) > 1 & length(ms2sub) <= 30){
  for(i in 1:length(ms2sub)){
    j <- order(intensitats)[i]
    
    raw_data <- readMSData(
      files = paste0("data/", polarity, "_DDA_mzmL/", ms2sub[[j]]@annotation), 
      mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = mz + 0.01 * c(-1, 1), 
                        rt = rt + 50 * c(-1, 1)
    )
    plot(chr, xlim = rt + 50 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    
    specplot(ms2sub[[j]],
             main = paste("id:", ms2sub[[j]]@id, " - ", 
                          ms2sub[[j]]@annotation))
    print(paste0(j, ": ", gsub(".*\\/", "", ms2sub[[j]]@annotation), " - ", ms2sub[[j]]@id, 
                 " - ", ms2sub[[j]]@rt))
  }
} else if(length(ms2sub) == 1){
  if(substr(ms2sub@annotation, 1, 1) == "x"){
    raw_data <- readMSData(files = paste0("data/", polarity, "_DDA_mzML/", ms2sub@annotation), 
                           mode = "onDisk")
  }
  chr <- chromatogram(raw_data, 
                      mz = mz + 0.01 * c(-1, 1), 
                      rt = rt + 20 * c(-1, 1)
  )
  plot(chr, xlim = rt + 20 * c(-1, 1))
  abline(v=ms2sub@rt)
  
  specplot(ms2sub,
           main = paste("id:", ms2sub@id, " - ", 
                        gsub("_DDA.*", "", gsub(".mzML", "", gsub(".*\\/", "", ms2sub@annotation)))),
           sub = "")
  print(paste0(gsub(".*\\/", "", ms2sub@annotation), " - ", ms2sub@id, " - ", ms2sub@rt))
}else {
  raw_data <- readMSData(files = ms2sub@annotation, mode = "onDisk")
  chr <- chromatogram(raw_data, 
                      mz = mz + 0.01 * c(-1, 1), 
                      rt = rt + 20 * c(-1, 1)
  )
  plot(chr, xlim = rt + 20 * c(-1, 1))
  abline(v=ms2sub@rt)
  specplot(ms2sub)
  print(ms2sub)
  tmp <- data.frame(ms2sub@spectrum)
  tmp <- tmp[tmp$X2 > tmp$X2[which.max(tmp$X2)]*0.01,  ]
}
