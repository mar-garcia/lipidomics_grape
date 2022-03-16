library(xcms)
library(MsCoreUtils)
load("data/RData.RData")


s <- "maturation" # specify "maturation" or "tissues"
p <- "POS" # specify "POS" or "NEG"
ft_mz <- 820.7418           

dda_xdata <- get(paste("dda_xdata", s, p, sep = "_"))
ft <- get(paste("ft", s, p, sep = "_"))

dt <- data.frame(t(featureValues(dda_xdata, value = "into")))
(ft_id <- rownames(featureDefinitions(dda_xdata, mz = ft_mz, ppm = 0.1)))
ft[ft_id,]
ft_id <- ft_id[2]

ydata <- filterFile(dda_xdata, which.max(dt[,ft_id]))
ydata <- ydata[msLevel(ydata) == 1]
idx <- which(rownames(ft) == ft_id)
chr <- chromatogram(ydata, mz = ft_mz + 0.01 * c(-1,1), 
                    ft$rtmed[idx] + 30 * c(-1, 1))
plot(chr)
myrt <- ft$rtmed[idx]
abline(v=myrt)

sps <- ydata[[closest(myrt, rtime(ydata))]]
sps <- cbind(mz(sps), intensity(sps))
plot(sps[,1], sps[,2], type = "h", xlim = ft_mz + 100 * c(-1,1))
idx <- which(sps[,2]/max(sps[,2]) > 0.1)
text(sps[idx,1], sps[idx,2], round(sps[idx,1], 4), cex = 0.8)


##############################
ft_chr1 <- featureChromatograms(dda_xdata, features = "FT182", 
                                expandRt = 15, filled = FALSE)
plot(ft_chr1)


##############################
library(CompoundDb)
cmps <- cbind(
  compounds(CompDb("data/CompDb.inhouse.0.sqlite"), "compound_id"), 
  compounds(CompDb("data/CompDb.inhouse.0.sqlite")))
cmps[grep(gsub("-", "|", "C0071-C0161-C0251-C0420-C1033-C1147-C2000-C2001-C2008"), cmps$compound_id),]
