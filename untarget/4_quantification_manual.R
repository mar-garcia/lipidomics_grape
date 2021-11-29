s <- "maturation"
cmp <- "PE(16:0/18:2)"

i <- which(res$compound == cmp)
fml <- unique(cmps$formula[cmps$compound == res$compound[i]])
cl <- unique(cmps$class[cmps$compound == res$compound[i]])
if(cl %in% c("PA", "mPA", "PE")){
  ad <- "[M-H]-"
  pol <- "NEG"
} else if(cl %in% c("PC", "MGDG", "DGDG")){
  ad <- "[M+CHO2]-"
  pol <- "NEG"
} else if(cl %in% c("DAG", "TAG")){
  ad <- "[M+NH4]+"
  pol <- "POS"
} else if(res$compound[i] %in% c("(Epi)catechin")){
  ad <- "[M-H]-"
  pol <- "NEG"
}
ft <- dbexp$FT[dbexp$compound == res$id[i] & dbexp$study == s]
xdata <- get(paste0("xdata_", s, "_", pol))
mz <- mass2mz(getMolecule(fml)$exactmass, ad)
ydata <- get(paste0("ydata_", s, "_", pol))
col <- get(paste0("col_", s, "_", pol))
if(length(ft) == 1){
  feat <- get(paste0("feat_", s, "_", pol))
  cp <- chromPeaks(xdata)[unlist(feat[ft, "peakidx"]),]
  rtmin <- min(cp[,"rtmin"])
  rtmax <- max(cp[,"rtmax"])
} else if(length(ft) == 0){
  cp <- as.data.frame(get(paste0("cp_", s, "_", pol)))
  cpx <- cp[unlist(matchWithPpm(mz, cp$mz, ppm = 10)),]
  cpx <- cpx[(cpx$rtmin < res$RT[i]*60) & (cpx$rtmax > res$RT[i]*60),]
  if(any(duplicated(cpx$sample))){
    print(paste("Ull! Check", res$compound[i], "in", s, "for duplicated cpx$sample!"))
  } else {
    rtmin <- min(cpx$rtmin)
    rtmax <- max(cpx$rtmax)
  }
} else {
  print(paste("ULL! Check", res$compound[i], "in", s))
}
i.cmp <- matrix(c(mz + 0.01 * c(-1, 1),
                  rtmin, rtmax), nrow = 1)
colnames(i.cmp) <- c("mzmin", "mzmax", "rtmin", "rtmax")
chr <- chromatogram(ydata, 
                    mz = i.cmp[,1:2], rt = i.cmp[,3:4] + 10 * c(-1,1))
plot(chr, col = col, xaxt = "n")
axis(1, at = seq(0, 60*30, 3), labels = sprintf("%.2f", seq(0, 30, 3/60))) 
abline(v=i.cmp[,3:4])