setwd("~/GitHub/lipidomics/untarget")
polarity <- "POS" # specify "POS" or "NEG"
load(paste0("data/RData/data_XCMS_", polarity, ".RData"))
data <- featureValues(xdata, method = "sum", value = "into")
features <- data.frame(featureDefinitions(xdata))

z.ft <- "FT033"
z.idx <- which(rownames(features) == z.ft)
# RT range
z.features <- features[(features$rtmed > (features$rtmed[z.idx] - 10)) & 
                         ( features$rtmed < (features$rtmed[z.idx] + 10)), ]
# intensity correlation
z.data <- data[rownames(data) %in% rownames(z.features), ]
z.data <- z.data[cor(t(z.data), z.data[z.ft,]) > 0.7,]
# peak-shape correlation
z.xdata <- filterFile(xdata, which(xdata$order == "x099"))
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
z.features <- z.features[z.cor > 0.7, ]

z.features



z.xdata %>%
  filterRt(rt = features$rtmed[z.idx] + 10 * c(-1,1)) %>%
  filterMz(mz = features$mzmed[z.idx]+1.003355*2 + 0.01 * c(-1, 1)) %>%
  plot(type = "XIC")
