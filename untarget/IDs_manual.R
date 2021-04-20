library(beeswarm)
load("ID_tissues.RData")
my.features$FGx[is.na(my.features$FGx)] <- "X"
data <- rbind(data_pos, data_neg)
class <- colnames(data)
class <- gsub("pt11_", "", class)
class <- gsub("_.*", "", class)
sp_ts <- grep("skin|pulp|seeds", colnames(data))
col_class_ts <- RColorBrewer::brewer.pal(8, "Set1")[6:8]#[c(1, 4, 5)]
names(col_class_ts) <- c("pulp", "seeds", "skin")
target <- readxl::read_xlsx("../target/data/data_tissues.xlsx", sheet = "raw_data")
colnames(target)[1] <- "ID"

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
  beeswarm(
    vals, col = col_class_ts[names(vals)], pch = 16,
    main = paste(rownames(my.features)[idx[i]], "-",
                 round(my.features$mzmed[rownames(my.features) == 
                                           rownames(my.features)[idx[i]]], 4)))
  bxplot(vals, probs = 0.5, col = "#00000060", add = TRUE)
}
# plot integrated areas
load("tissues/data/RData/data_XCMS_POS.RData")
chrs <- featureChromatograms(xdata, features = c("FT0116", "FT0118"))
sample_colors <- col_type[xdata$tissue]
sample_colors[!names(sample_colors) %in% c("seeds", "pulp", "skin")] <- "#999999"
plot(chrs, peakBg = paste0(sample_colors[chromPeaks(chrs)[, "sample"]], "60"))

# overplot EICs
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

# Correlation with target data
tg <- target[target$ID == "18:1(d7) Lyso PC",grep("pt11", colnames(target))]
ut <- data["P0116", grep("pt11", colnames(data))]
plot(t(tg), ut, col = col_type[class[!grepl("QC|xx00", class)]], pch = 16,
     main = paste("corr", round(cor(t(tg), ut), 3)), xlab = "target", ylab = "untarget")

#############################################################################

x.feat <- my.features[idx, ]
x.feat <- x.feat[order(x.feat$mzmed),]
x.feat[, c("mzmed", "rtmed", "FG", "MS2", "isotopes", "FGx", "annotation")]
feat[feat$X %in% rownames(x.feat),]




x.feat <- my.features[my.features$FGx != "X", ]
x.feat <- x.feat[order(x.feat$mean_max2, decreasing = T), ]
head(x.feat)
