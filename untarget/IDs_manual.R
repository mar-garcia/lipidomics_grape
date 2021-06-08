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
x <- "P0001"
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
chrs <- featureChromatograms(
  xdata, 
  features = gsub("P", "FT", rownames(my.features)[idx]))
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
legend("topright", pch = 16, col = seq(length(idx)), legend = round(my.features$mzmed[idx],4))


# Correlation with target data ----
target$ID[grep("PA\\(16:0", target$ID)]
tg <- target[target$ID == "TAG54:7_FA18:2",grep("pt11", colnames(target))]
ut <- data[rownames(my.features)[idx][1], grep("pt11", colnames(data))]
plot(t(tg), ut, col = col_class_ts[class[!grepl("QC|xx00", class)]], pch = 16,
     main = paste("corr", round(cor(t(tg), ut), 3)), xlab = "target", ylab = "untarget")
my.features[idx[1],"mzmed"]
my.features[idx[1],"rtmed"]/60

# MS2 ---------------
mz <- 984.8944 
rt <- 22.41*60
#ms2sub <- getSpectrum(ms2list, "precursor", mz, mz.tol = 0.001)
ms2sub <- filterPrecursorMz(ms2, mz + 0.01 * c(-1, 1))
#ms2sub <- getSpectrum(ms2sub, "rt", rt, rt.tol = 10)
ms2sub <- filterRt(ms2sub, rt + 10 * c(-1, 1))
if(length(ms2sub) > 1){
  intensitats <- c()
  for(i in seq(length(ms2sub))){
    idx <- c(which(unlist(mz(ms2sub[i])) > 283.2 & unlist(mz(ms2sub[i])) < 283.9),
             which(unlist(mz(ms2sub[i])) > 341.1 & unlist(mz(ms2sub[i])) < 341.3)
    )
    int.noise <- max(unlist(intensity(ms2sub[i]))[idx])
    int.good <- max(unlist(intensity(ms2sub[i]))[-idx])
    intensitats <- c(intensitats, int.good / int.noise)
  }
}
rm(tmp)
length(intensitats)
#tmp <- order(intensitats)[(length(intensitats)-10):length(intensitats)]
dev.off()
if(exists("tmp")){
  i.seq <- tmp
} else {
  i.seq <- seq(length(ms2sub))
}
par(mfrow = c(1, 2))
for(i in i.seq){
  if(exists("tmp")){
    j <- i
  } else {
    j <- order(intensitats)[i]
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
}




# Sirius ----
load("tissues/data/RData/data_XCMS_POS.RData")
y.xdata <- filterFile(xdata, which.max(data[rownames(my.features)[idx[1]],]))
mz <- 869.7390 
rt <- 20.70*60
xdata <- readMSData(
  files = fileNames(y.xdata),
  mode = "onDisk")
chr <- chromatogram(xdata, mz = mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
rt2 <-  1242.60                     
plot(chr, xlim = c(rt2 - 50, rt2 + 50))
abline(v = rt2)
sps <- xdata[[closest(rt2, rtime(xdata)#, duplicates = "closest"
)]]
sps <- as.data.frame(sps)
plot(sps$mz, sps$i, type = "h", xlim = c(mz - 10, mz + 10))
text(sps$mz, sps$i, round(sps$mz, 4), cex = 0.8)
sps[unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)):
      (unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10))+4),]
tmp <- sps[((unlist(CompoundDb::matchWithPpm(mz, sps$mz, ppm = 10)))-4):
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
