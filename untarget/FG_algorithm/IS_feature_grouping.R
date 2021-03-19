setwd("~/GitHub/lipidomics/untarget/tissues")
library(xcms)
library(CluMSID)

z.polarities <- c("POS", "NEG")
for(z in seq(length(z.polarities))){
  load(paste0("data/RData/IS_data_XCMS_", z.polarities[z], ".RData"))
  load(paste0("data/RData/MS2_library_", z.polarities[z], ".RData"))
  data <- featureValues(xdata, method = "sum", value = "into")
  data[is.na(data)] <- 0
  features <- data.frame(featureDefinitions(xdata))
  features$polarity <- z.polarities[z]
  features$mean_solv <- rowMeans(data[,xdata$tissue == "solv"])
  features$mean_blank <- rowMeans(data[,xdata$tissue == "blank"])
  features$mean_STDmix <- rowMeans(data[,xdata$tissue == "STDmix"])
  features$mean_max <- apply(
    features[,c("mean_solv", "mean_blank", "mean_STDmix")], 1, max)
  features <- features[order(features$mean_max, decreasing = T), ]
  features <- features[features$npeaks == 2 & features$blank == 2 & 
                         features$solv == 0 & features$STDmix == 0, ]
  features$FG <- NA
  y.n <- 0
  for(y in seq(nrow(features))){
    if(is.na(features$FG[y])){
      y.n <- y.n + 1 
      y.features <- features[is.na(features$FG), ]
      y.ft <- which(rownames(y.features) == rownames(features)[y])
      # RT range
      y.features <- y.features[(y.features$rtmed > (y.features$rtmed[y.ft] - 10)) & 
                                 ( y.features$rtmed < (y.features$rtmed[y.ft] + 10)), ]
      if(nrow(y.features) > 1){
        # intensity correlation
        y.data <- data[rownames(data) %in% rownames(y.features), ]
        y.cor <- cor(t(y.data), y.data[rownames(features)[y],])
        y.features <- merge(y.features,y.cor, by = "row.names")
        colnames(y.features)[ncol(y.features)] <- "corr_int"
        rownames(y.features) <- y.features$Row.names
        y.features <- y.features[,-1]
        y.features <- y.features[y.features$corr_int > 0.7, ]
        
        # peak-shape correlation
        y.xdata <- filterFile(xdata, #which(xdata$order == "x099")
                              which.max(data[rownames(features)[y],]))
        y.chr <- chromatogram(y.xdata,
                              mz = features$mzmed[y] + 0.01 * c(-1, 1),
                              rt = features$rtmed[y] + 10 * c(-1,1),
                              aggregationFun = "max")
        y.cor <- c()
        for(i in seq(nrow(y.features))){
          y.chr2 <- chromatogram(y.xdata,
                                 mz = y.features$mzmed[i] + 0.01 * c(-1, 1),
                                 rt = y.features$rtmed[i] + 10 * c(-1,1),
                                 aggregationFun = "max")
          y.cor <- c(y.cor, correlate(y.chr[[1]], y.chr2[[1]]))
          names(y.cor)[length(y.cor)] <- rownames(y.features)[i]
        }
        y.cor[is.na(y.cor)] <- 0
        y.features <- merge(y.features,y.cor, by = "row.names")
        colnames(y.features)[ncol(y.features)] <- "corr_ps"
        rownames(y.features) <- y.features$Row.names
        y.features <- y.features[,-1]
        
        # MS2 associations
        y.features$MS2 <- NA
        for(i in seq(nrow(y.features))){
          y.ms2 <- findFragment(ms2list, y.features$mzmed[i])
          y.ms2 <- getSpectrum(y.ms2, "rt", y.features$rtmed[i], rt.tol = 10)
          for(j in seq(nrow(y.features))){
            if(length(y.ms2) == 1){
              tmp <- abs(y.ms2@precursor - y.features$mzmed[j]) < 0.01
            } else if(length(y.ms2) > 1){
              tmp <- length(getSpectrum(y.ms2, "precursor", 
                                        y.features$mzmed[j], mz.tol = 0.01)) > 0
            } else if(length(y.ms2) == 0){
              tmp <- FALSE
            }
            if(is.na(y.features$MS2[i])){
              y.features$MS2[i] <- tmp
            } else {
              y.features$MS2[i] <- paste0(y.features$MS2[i], "; ", tmp)
            }
            
          }
        }
        y.features$MS2 <- grepl("TRUE", y.features$MS2)
        
        y.features <- y.features[
          y.features$corr_ps > 0.7 |
            (y.features$corr_ps > 0.5 & y.features$MS2), ]  
      }
      features$FG[rownames(features) %in% rownames(y.features)] <- y.n
    }
  }
  
  
  
  features$FG <- paste0(substr(z.polarities[z], 1, 1), 
                        formatC(features$FG, width = nchar(nrow(features)), 
                                flag = "0"))
  if(z == 1){
    featuresz <- features
  } else if(z == 2){
    features <- rbind(featuresz, features)
  }
}
rm(z, xdata, ms2list, data, featuresz, z.polarities,
   y, y.n, y.ft, y.features, y.data, y.xdata, y.chr, y.chr2, y.cor, i, 
   j, y.ms2, tmp)


features <- features[order(features$mean_max, decreasing = T), ]
features$FGx <- NA
y.n <- 0
y <- 1
if(is.na(features$FGx[y])){
  y.features <- features[is.na(features$FGx), ]
  y.ft <- which(rownames(y.features) == rownames(features)[y])
  # Within the same polarity, take just only the features aggruped togheter "y", 
  # whereas keep all the other features of the other polarity
  y.features <- y.features[(y.features$FG == y.features$FG[y.ft]) | 
                             (y.features$polarity != y.features$polarity[y.ft]),]
  # RT range
  y.features <- y.features[(y.features$rtmed > (y.features$rtmed[y.ft] - 10)) & 
                             ( y.features$rtmed < (y.features$rtmed[y.ft] + 10)), ]
}
