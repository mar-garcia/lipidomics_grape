setwd("~/GitHub/lipidomics/untarget")
polarity <- "POS" # specify "POS" or "NEG"
load(paste0("data/RData/data_XCMS_", polarity, ".RData"))
data <- featureValues(xdata, method = "sum", value = "into")


ft <- "FT14"
which.max(data[ft,])
