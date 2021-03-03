library(xcms)

polarity <- "POS" # specify "POS" or "NEG"
load(paste0("data/RData/data_XCMS_", polarity, ".RData"))

metadata <- pData(xdata)
data <- featureValues(xdata, method = "sum", value = "into")
data[is.na(data)] <- 0
data <- data.frame(t(data))

data <- data[metadata$type == "study", ]
metadata <- metadata[metadata$type == "study", ]
metadata$time <- as.numeric(metadata$time)


fun.lm <- function(x) {
  analyte <- data[,x]
  lmreg <- lm(analyte ~ metadata$time,
              data = data)
  res <- coef(summary(lmreg))[-1, ]
  p.values <- data.frame(c(res[[1]],
                           res[[4]]))
  names(p.values) <- x
  rownames(p.values) <- c(paste0("coef_", rownames(res)),
                          paste0("p-value_", rownames(res)))
  return(p.values)
}
lm.tests <- do.call(cbind.data.frame,
                      lapply(colnames(data[,]),
                             function(x) fun.lm(x)))
lm.tests <- data.frame(t(lm.tests))
lm.tests$padj <- p.adjust(lm.tests$p.value_, "bonferroni")
sum(lm.tests$padj < 0.05)
i <- 2
plot(metadata$time, data[,which(lm.tests$padj < 0.05)[i]], 
     xlab = "time", ylab = colnames(data)[which(lm.tests$padj < 0.05)[i]])
