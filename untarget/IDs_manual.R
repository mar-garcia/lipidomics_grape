library(beeswarm)
library(xcms)
library(CompoundDb)
library(CluMSID)
library(MsCoreUtils)
library(Spectra)
load("ID_tissues.RData")
my.features$FGx[is.na(my.features$FGx)] <- "X"
data <- rbind(data_pos, data_neg)
class <- colnames(data)
class <- gsub("pt11_", "", class)
class <- gsub("_.*", "", class)
col_type <- RColorBrewer::brewer.pal(9, name = "Set1")
if(study == "maturation"){
  names(col_type) <- c("xQC",      # red
                       "slv",      # blue
                       "study",    # green
                       "QCrw_05uL",# viola
                       "QCdl_05uL"#,# orange
                       #"slveq",   # yellow
                       #"QCeq",    # brown
                       #"STDmix"   # pink
  )
  col_type <- c(col_type[1:3],
                rep(col_type[4], 3),
                rep(col_type[5], 3))
  names(col_type) <- c("xQC", "slv", "study", 
                       "QCrw_02uL", "QCrw_05uL", "QCrw_10uL", 
                       "QCdl_02uL", "QCdl_05uL", "QCdl_10uL")
} else if(study == "tissues"){
  names(col_type) <- c("QC",      # red
                       "solv",      # blue
                       "entire",    # green
                       "QCrw_05uL",# viola
                       "QCdl_05uL",# orange
                       "pulp",   # yellow
                       "seeds",    # brown
                       "skin"   # pink
  )
  col_type <- c(col_type[1:2], col_type[c(3,6:8)],
                rep(col_type[4], 3),
                rep(col_type[5], 3))
  names(col_type) <- c("QC", "solv", "entire", "pulp", "seeds", "skin", 
                       "QCrw_02uL", "QCrw_05uL", "QCrw_10uL", 
                       "QCdl_02uL", "QCdl_05uL", "QCdl_10uL")
}
sp_ts <- grep("skin|pulp|seeds", colnames(data))
col_class_ts <- RColorBrewer::brewer.pal(8, "Set1")[c(5, 4, 3)]
names(col_class_ts) <- c("seeds", "skin", "pulp")
target <- readxl::read_xlsx("../target/data/data_tissues.xlsx", sheet = "raw_data")
colnames(target)[1] <- "ID"


# Compound ----
x <- "P0002"
idx <- which(my.features$FGx == x)
if(length(idx == 2)){
  par(mfrow = c(1, 2))
} else{
  par(mfrow = c(3, 2))
}
for(i in seq(length(idx))){
  vals <- split(data[rownames(my.features)[idx[i]], sp_ts],
                f = class[sp_ts])
  vals <- vals[c("skin", "pulp", "seeds")]
  beeswarm(
    vals, col = col_class_ts[names(vals)], pch = 16,
    main = paste(rownames(my.features)[idx[i]], "-",
                 round(my.features$mzmed[rownames(my.features) == 
                                           rownames(my.features)[idx[i]]], 4)))
  bxplot(vals, probs = 0.5, col = "#00000060", add = TRUE)
}

# plot integrated areas ------
load("tissues/data/RData/data_XCMS_POS.RData")
chrs <- featureChromatograms(xdata, features = c("FT0116", "FT0118"))
sample_colors <- col_type[xdata$tissue]
sample_colors[!names(sample_colors) %in% c("seeds", "pulp", "skin")] <- "#999999"
plot(chrs, peakBg = paste0(sample_colors[chromPeaks(chrs)[, "sample"]], "60"))

load("tissues/data/RData/data_XCMS_3ts_POS.RData")
ft <- data.frame(featureDefinitions(xdata))
idx2 <- c()
for(i in seq(length(idx2))){
  idx2 <- c(idx2, rownames(ft)[unlist(matchWithPpm(my.features$mzmed[
    rownames(my.features) == rownames(my.features)[idx[i]]], ft$mzmed, ppm = 10))])
}
chrs <- featureChromatograms(xdata, features = idx2)
sample_colors <- col_class_ts[xdata$tissue]
plot(chrs, peakBg = paste0(sample_colors[chromPeaks(chrs)[, "sample"]], "60"))

# overplot EICs ----
load("tissues/data/RData/data_XCMS_POS.RData")
par(mfrow=c(1, 2))
i <- 1
y.xdata <- filterFile(xdata, which.max(data[rownames(my.features)[idx[i]],]))
y.chr <- chromatogram(y.xdata,
                      mz = my.features$mzmed[rownames(my.features) == 
                                               rownames(my.features)[idx[i]]] + 0.01 * c(-1, 1),
                      rt = my.features$rtmed[rownames(my.features) == 
                                               rownames(my.features)[idx[i]]] + 10 * c(-1,1),
                      aggregationFun = "max", include = "none")
plot(y.chr)
for(i in 2:length(idx)){
  y.chr2 <- chromatogram(y.xdata,
                         mz = my.features$mzmed[rownames(my.features) == 
                                                  rownames(my.features)[idx[i]]] + 0.01 * c(-1, 1),
                         rt = my.features$rtmed[rownames(my.features) == 
                                                  rownames(my.features)[idx[i]]] + 10 * c(-1,1),
                         aggregationFun = "max", include = "none")
  points(rtime(y.chr2[[1]]), intensity(y.chr2[[1]]), col = i, type = "l")
}

i <- 1
y.xdata <- filterFile(xdata, which.max(data[rownames(my.features)[idx[i]],]))
plot(rtime(y.chr[[1]]), intensity(y.chr[[1]])/max(intensity(y.chr[[1]]), na.rm = T), col = i, type = "l")
for(i in 2:length(idx)){
  y.chr2 <- chromatogram(
    y.xdata,
    mz = my.features$mzmed[rownames(my.features) == 
                             rownames(my.features)[idx[i]]] + 0.01 * c(-1, 1),
    rt = my.features$rtmed[rownames(my.features) == 
                             rownames(my.features)[idx[i]]] + 10 * c(-1,1),
    aggregationFun = "max", include = "none")
  points(rtime(y.chr2[[1]]), intensity(y.chr2[[1]])/max(intensity(y.chr2[[1]]), na.rm = T), 
         col = i, type = "l")
}
legend("topright", pch = 16, col = seq(length(idx)), legend = round(my.features$mzmed[rownames(my.features) == 
                                                                                        rownames(my.features)[idx]],4))

# Correlation with target data ----
tg <- target[target$ID == "18:1(d7) Lyso PC",grep("pt11", colnames(target))]
ut <- data["P0116", grep("pt11", colnames(data))]
plot(t(tg), ut, col = col_class_ts[class[!grepl("QC|xx00", class)]], pch = 16,
     main = paste("corr", round(cor(t(tg), ut), 3)), xlab = "target", ylab = "untarget")

# MS2
mz <- 573.3903
rt <- 8.27*60
ms2sub <- getSpectrum(ms2list, "precursor", mz, mz.tol = 0.001)
ms2sub <- getSpectrum(ms2sub, "rt", rt, rt.tol = 10)
if(length(ms2sub) > 1){
  intensitats <- c()
  for(i in seq(ms2sub)){
    idx <- which(accessSpectrum(ms2sub[[i]])[,1] > 283.2 & accessSpectrum(ms2sub[[i]])[,1] < 283.9)
    if(length(idx) == 0){
      idx <- substring(gsub(".*\\.","", accessSpectrum(ms2sub[[i]])[,1]), 1, 1)>1 
    }
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
    
    if(grepl("lipidgrape_tissues", ms2sub[[j]]@annotation)){
      study <- "tissues"
    } else{
      study <- "maturation"
    }
    
    raw_data <- readMSData(
      files = paste0(study, "/data/", gsub(".*_", "", gsub("_DDA.mzML", "", ms2sub[[j]]@annotation)), "_DDA_mzmL/", ms2sub[[j]]@annotation), 
      mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = mz + 0.01 * c(-1, 1), 
                        rt = rt + 20 * c(-1, 1)
    )
    plot(chr, xlim = rt + 20 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    
    specplot(ms2sub[[j]], main = ms2sub[[j]]@id)
    print(paste0(j, ": ", gsub(".*\\/", "", ms2sub[[j]]@annotation), " - ", ms2sub[[j]]@id, 
                 " - ", ms2sub[[j]]@rt))
  }
} else if(length(ms2sub) > 1 & length(ms2sub) <= 30){
  for(i in 1:length(ms2sub)){
    j <- order(intensitats)[i]
    
    if(grepl("lipidgrape_tissues", ms2sub[[j]]@annotation)){
      study <- "tissues"
    } else{
      study <- "maturation"
    }
    
    
    raw_data <- readMSData(
      files = paste0(study, "/data/", gsub(".*_", "", gsub("_DDA.mzML", "", ms2sub[[j]]@annotation)), "_DDA_mzmL/", ms2sub[[j]]@annotation), 
      mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = mz + 0.01 * c(-1, 1), 
                        rt = rt + 50 * c(-1, 1)
    )
    plot(chr, xlim = rt + 50 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    
    specplot(ms2sub[[j]], main = ms2sub[[j]]@id)
    print(paste0(j, ": ", gsub(".*\\/", "", ms2sub[[j]]@annotation), " - ", ms2sub[[j]]@id, 
                 " - ", ms2sub[[j]]@rt))
  }
} else if(length(ms2sub) == 1){
  if(grepl("lipidgrape_tissues", ms2sub@annotation)){
    study <- "tissues"
  } else{
    study <- "maturation"
  }
  if(substr(ms2sub@annotation, 1, 1) == "x"){
    raw_data <- readMSData(files = paste0(study, "/data/", polarity, "_DDA_mzML/", ms2sub@annotation), 
                           mode = "onDisk")
  }
  chr <- chromatogram(raw_data, 
                      mz = mz + 0.01 * c(-1, 1), 
                      rt = rt + 20 * c(-1, 1)
  )
  plot(chr, xlim = rt + 20 * c(-1, 1))
  abline(v=ms2sub@rt)
  
  specplot(ms2sub, main = ms2sub@id)
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

# Sirius ----
xdata <- readMSData(
  files = fileNames(y.xdata),
  mode = "onDisk")
chr <- chromatogram(xdata, mz = mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
rt2 <- 499.018                    
plot(chr, xlim = c(rt2 - 50, rt2 + 50))
abline(v = rt2)
sps <- xdata[[closest(rt2, rtime(xdata)#, duplicates = "closest"
)]]
sps <- as.data.frame(sps)
plot(sps$mz, sps$i, type = "h", xlim = c(mz - 10, mz + 10))
text(sps$mz, sps$i, round(sps$mz, 4), cex = 0.8)
sps[unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)):
      (unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10))+4),]
tmp <- sps[unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)):
             (unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10))+4),]
tmp  <- tmp[order(tmp$i, decreasing = T), ]
tmp[1:4,]
write.table(tmp[1:3,], 
            paste0("identifications/sirius/", round(mz), "_", 
                   round(rt2), "_FS.txt"), row.names = F, col.names = F)

ms2sub <- getSpectrum(ms2list, "annotation", "x010_lipidgrape_tissues_pt11_skin_rep3_NEG_DDA.mzML")
ms2sub <- getSpectrum(ms2sub, "precursor", mz, mz.tol = 0.001)
ms2sub <- getSpectrum(ms2sub, "rt", 497, rt.tol = 1)
write.table(ms2sub@spectrum, 
            paste0("identifications/sirius/", round(mz), "_", 
                   round(rt2), "_MS2.txt"), 
            row.names = F, col.names = F)

# Spectra ----
sps <- Spectra("tissues/data/NEG_DDA_mzML/x010_lipidgrape_tissues_pt11_skin_rep3_NEG_DDA.mzML")
sps <- filterPrecursorMz(sps, mz + 0.01 * c(-1,1))
filterRt(sps, 497 + 1 * c(-1,1))

# Others ----

x.feat <- my.features[idx, ]
x.feat <- x.feat[order(x.feat$mzmed),]
x.feat[, c("mzmed", "rtmed", "FG", "MS2", "isotopes", "FGx", "annotation")]
feat[feat$X %in% rownames(x.feat),]




x.feat <- my.features[my.features$FGx != "X", ]
x.feat <- x.feat[order(x.feat$mean_max2, decreasing = T), ]
head(x.feat)
