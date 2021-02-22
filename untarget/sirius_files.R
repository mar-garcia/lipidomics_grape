library(CluMSID)
library(MsCoreUtils)
setwd("~/GitHub/lipidomics/untarget")
polarity <- "POS" # specify "POS" or "NEG"
load(paste0("data/RData/data_XCMS_", polarity, ".RData"))
load(paste0("data/RData/MS2_library_", polarity, ".RData"))

data <- featureValues(xdata, method = "sum", value = "into")
features <- data.frame(featureDefinitions(xdata))

ft <- "FT04"
mz <- features$mzmed[rownames(features)==ft]
rt <- features$rtmed[rownames(features)==ft]

which.max(data[ft,])
xdata <- readMSData(
  files = paste0("data/", polarity, "_FS_fixed/", names(which.max(data[ft,]))),
  mode = "onDisk")
chr <- chromatogram(xdata, mz = mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
rt2 <- 1070.680
plot(chr, xlim = c(rt2 - 50, rt2 + 50))
abline(v = rt2)
sps <- xdata[[closest(rt2, rtime(xdata)#, duplicates = "closest"
)]]
sps <- as.data.frame(sps)
plot(sps$mz, sps$i, type = "h", xlim = c(mz - 10, mz + 10))
text(sps$mz, sps$i, round(sps$mz, 4), cex = 0.8)
sps[unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)):
    (unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10))+9),]
tmp <- sps[unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)):
             (unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10))+9),]
tmp  <- tmp[order(tmp$i, decreasing = T), ]
tmp[1:8,]
write.table(tmp[1:8,], paste0("data/sirius/", round(mz), "_", round(rt2), "_FS.txt"), row.names = F, col.names = F)


##################################################################

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
    
    raw_data <- readMSData(files = ms2sub[[j]]@annotation, mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = mz + 0.01 * c(-1, 1), 
                        rt = rt + 50 * c(-1, 1)
    )
    plot(chr, xlim = rt + 50 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    
    specplot(ms2sub[[j]],
             main = paste("id:", ms2sub[[j]]@id, " - ", 
                          gsub(paste0("_DDA_", polarity, ".mzML"), "", 
                               gsub("_DDA.*", "", gsub(".mzML", "", gsub(".*\\/", "", ms2sub[[j]]@annotation))))))
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
i <- 5
write.table(ms2sub[[i]]@spectrum, paste0("data/sirius/", ms2sub[[j]]@id, "_MS2.txt"), row.names = F, col.names = F)
