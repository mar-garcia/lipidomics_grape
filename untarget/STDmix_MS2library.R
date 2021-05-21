library(readxl)
library(Spectra)
library(pheatmap)
library(MsCoreUtils)
library(CompoundDb)
library(OrgMassSpecR)
library(Rdisop)
library(igraph)
library(xcms)

norm_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}
plotms2 <- function(x){
  mz <- unlist(mz(x))
  int <- unlist(intensity(x))
  i.thr <- 10
  plot(mz, int, type = "h", xlab = "mz", ylab = "relative intensity",
       main = paste(sprintf("%.4f", x@backend@spectraData$precursorMz), 
                    "@", sprintf("%.2f", x@backend@spectraData$rtime/60)))
  idx <- which(int > i.thr)
  text(mz[idx], int[idx], sprintf("%.4f", mz[idx]), cex = 0.8)
}

plotms2.chr <- function(x){
  xdata <- readMSData(x@backend@spectraData$dataOrigin, mode = "onDisk")
  chr <- chromatogram(
    xdata, mz = x@backend@spectraData$precursorMz + 0.01 * c(-1, 1),
    rt = x@backend@spectraData$rtime + 30 * c(-1, 1)
  )
  plot(chr, xaxt="n",
       main = gsub("_DDA.mzML", "", basename(x@backend@spectraData$dataOrigin)))
  axis(1, at = seq(0, 60*30, 6), labels = sprintf("%.2f", seq(0, 30, 6/60))) 
  points(x@backend@spectraData$rtime, 
         intensity(chr[[1]])[closest(x@backend@spectraData$rtime, 
                                     rtime(chr[[1]]))], pch = 8) 
}

#####################################################


cmps <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/db_compounds.xlsx")
cmps[, 5:ncol(cmps)] <- cmps[,5:ncol(cmps)] * 60
inj <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/standards_injections.xlsx")
inj <- inj[!is.na(inj$ID_cmp),]

# get RT only from lipidomics chromatography
cmps$RT <- rowMeans(cmps[,c("lipidgrape_min", "lipidgrape_max")])
cmps <- subset(cmps, select = c("ID_cmp", "compound", "formula", "RT"))

# select Lyso-PC
cmpsx <- cmps[grep("TAG", cmps$compound),]

p <- "[M+NH4]+"
d <- 0
i.pol <- "POS"
j <- 9
i.inj <- inj[grep(paste(cmpsx$ID_cmp, collapse = "|"), inj$ID_cmp),]
i.inj <- i.inj[i.inj$polarity == i.pol, ]
i.inj <- i.inj[grep(cmpsx$ID_cmp[j], i.inj$ID_cmp),]
i.fml <- cmpsx$formula[j]
rm(sps.ms2)
for(i in seq(nrow(i.inj))){
  i.fl <- Spectra(paste0("C:/Users/garciaalom/Documents/GitHub/MS2_library/mzML/",
                         i.inj$path[i], i.inj$filename[i], ".mzML"), 
                  backend = MsBackendDataFrame())
  i.fl <- i.fl[msLevel(i.fl) > 1]
  if(grepl("D", i.fml)){
    i.mz <- unlist(mass2mz(getMolecule(i.fml)$exactmass, p))
  } else {
    i.mz <- unlist(mass2mz(MonoisotopicMass(formula = ListFormula(i.fml)), p))
  }
  i.mz <- i.mz + d
  i.fl <- filterPrecursorMz(i.fl, i.mz + 0.01 * c(-1, 1))
  i.rt <- cmpsx$RT[j]
  i.fl <- filterRt(i.fl, i.rt + (10 * c(-1, 1)))
  if(length(i.fl) > 0){
    i.fl$name <- cmpsx$ID_cmp[j]
    if(exists("sps.ms2")){
      sps.ms2 <- c(sps.ms2, i.fl)
    } else{
      sps.ms2 <- i.fl
    } # end if(exists("sps.ms2")
  }
} # end file "i"

length(sps.ms2)
sps.ms2 <- addProcessing(sps.ms2, norm_int)
res <- Spectra::compareSpectra(sps.ms2, ppm = 20)
#colnames(res) <- rownames(res) <- sps.ms2$name
hm <- pheatmap(res)
par(mfrow = c(3, 2))
for(i in hm$tree_row$order){plotms2(sps.ms2[i])}
idx <- seq(9)
#for(i in idx){
#  k <- hm$tree_row$order[i]
#  plotms2(sps.ms2[k])
#}
sps <- Spectra::combineSpectra(sps.ms2[hm$tree_row$order[idx]],
                      #f = pe_spectra$name,
                      #p = ms2_spectra$CLUSTER_ID,
                      intensityFun = base::mean,
                      mzFun = base::mean,
                      tolerance = 0.005,
                      #ppm = 0,
                      minProp = 0.5,
                      peaks = "intersect",
                      weighted = TRUE)
par(mfrow = c(1, 1))
plotms2(sps)
tmp <- c()
for(i in idx){
  k <- hm$tree_row$order[i]
  tmp <- c(tmp, paste(i.inj$ID_file[i.inj$filename == gsub(".mzML", "", basename(sps.ms2[k]@backend@spectraData$dataOrigin))],
                      basename(sps.ms2[k]@backend@spectraData$dataOrigin), 
                      sps.ms2[k]@backend@spectraData$acquisitionNum))
}
cmpsx$ID_cmp[j]
cmpsx$compound[j]
tmp[order(tmp)]
write.csv(tmp[order(tmp)], "x.csv")

###########################################################################

ms2lib <- read_xlsx("MS2_library.xlsx")
#i.pol <- "POS"
#i.inj <- inj[inj$polarity == i.pol, ]
#cmpsx <- cmps[grep("Lyso PC|Lyso_PC", cmps$compound),]
i.inj <- inj[inj$ID_file %in% ms2lib$ID_file, ]


cmps <- cmps[grep("PE\\(", cmps$compound), ]
cmps <- cmps[cmps$compound != "PE(16:0/18:3)",]
ms2lib <- ms2lib[ms2lib$ID_cmp %in% cmps$ID_cmp, ]
ms2lib <- ms2lib[grep("\\-", ms2lib$adduct),]

tmp <- unique(paste(ms2lib$ID_cmp, ms2lib$adduct))
rm(sps.ms2)
for(i in seq(length(tmp))){
  i.ms2lib <- ms2lib[paste(ms2lib$ID_cmp, ms2lib$adduct) == tmp[i], ]
  i.ms2lib <- i.ms2lib[1,]
  i.fls <- unique(i.ms2lib$ID_file)
  for(j in seq(length(i.fls))){
    idx <- which(i.inj$ID_file == i.fls[j])
    sps <- Spectra(paste0("C:/Users/garciaalom/Documents/GitHub/MS2_library/mzML/",
                          i.inj$path[idx], i.inj$filename[idx], ".mzML"), 
                   backend = MsBackendDataFrame())
    sps <- filterAcquisitionNum(
      sps, 
      as.integer(i.ms2lib$acquisitionNum[i.ms2lib$ID_file == i.fls[j]]))
    if(exists("spsx")){
      spsx <- c(spsx, sps)
    } else{
      spsx <- sps
    }
  }
  spsx <- addProcessing(spsx, norm_int)
  sps <- Spectra::combineSpectra(spsx,
                        #f = pe_spectra$name,
                        #p = ms2_spectra$CLUSTER_ID,
                        intensityFun = base::mean,
                        mzFun = base::mean,
                        tolerance = 0.005,
                        #ppm = 0,
                        minProp = 0.5,
                        peaks = "intersect",
                        weighted = TRUE)
  rm(spsx)
  sps$C <- i.ms2lib$ID_cmp[1]
  sps$name <- cmps$compound[cmps$ID_cmp == i.ms2lib$ID_cmp[1]]
  sps$adduct <- i.ms2lib$adduct[1]
  if(exists("sps.ms2")){
    sps.ms2 <- c(sps.ms2, sps)
  } else{
    sps.ms2 <- sps
  } # end if(exists("sps.ms2")
}


sps.ms2 <- replaceIntensitiesBelow(sps.ms2, threshold = 10, value = 0)
sps.ms2 <- filterIntensity(sps.ms2, intensity = c(0.1, 100))
res <- Spectra::compareSpectra(sps.ms2, ppm = 20)
colnames(res) <- rownames(res) <- paste(sps.ms2$name, sps.ms2$adduct)
hm <- pheatmap(res)




sps.ms2_nl <- applyProcessing(sps.ms2)
mz(sps.ms2_nl@backend) <- mz(sps.ms2_nl) - precursorMz(sps.ms2_nl)
res <- Spectra::compareSpectra(sps.ms2_nl, tolerance = 0.005)
colnames(res) <- rownames(res) <- paste(sps.ms2$name, sps.ms2$adduct)
hm <- pheatmap(res)






frags <- data.frame(table(round(unlist(mz(sps.ms2)), 3)))
frags <- frags[frags$Freq > 1, ]
frags$Var1 <- as.numeric(as.character(frags$Var1))
dt_f <- data.frame(matrix(ncol = nrow(frags), nrow = length(sps.ms2)))
colnames(dt_f) <- paste0("i_", frags$Var1)
rownames(dt_f) <- paste(sps.ms2$name, sps.ms2$adduct)
for(i in seq(nrow(frags))){
  dt_f[,i] <- ifelse(containsMz(sps.ms2, frags$Var1[i], 
                                tolerance = 0.01), 1, 0)
}
idx <- c()
for(i in seq(ncol(dt_f)-1)){
  if(identical(dt_f[,i], dt_f[,i+1]) &
     length(unlist(matchWithPpm(frags$Var1[i], frags$Var1[i+1], ppm = 10))) > 0){
    idx <- c(idx, i)
  }
}
dt_f <- dt_f[,-idx]

nl <- data.frame(table(round(unlist(mz(sps.ms2_nl)), 3)))
nl <- nl[nl$Freq > 1, ]
nl$Var1 <- as.numeric(as.character(nl$Var1))
dt_nl <- data.frame(matrix(ncol = nrow(nl), nrow = length(sps.ms2)))
colnames(dt_nl) <- paste0("l_", -as.numeric(nl$Var1))
rownames(dt_nl) <- paste(sps.ms2$name, sps.ms2$adduct)
for(i in seq(nrow(nl))){
  dt_nl[,i] <- ifelse(containsMz(sps.ms2_nl, nl$Var1[i], 
                                 tolerance = 0.01), 1, 0)
}
idx <- c()
for(i in seq(ncol(dt_nl)-1)){
  if(identical(dt_nl[,i], dt_nl[,i+1]) &
     abs(nl$Var1[i] - nl$Var1[i+1]) < 0.01){
    idx <- c(idx, i)
  }
}
dt_nl <- dt_nl[,-idx]

dt <- cbind(dt_f, dt_nl)
dt <- dt[,dt["PA(16:0/18:2) [M-H]-",] > 0]
dt <- dt[rowSums(dt) > 0,]
dt <- dt[, colSums(dt) > 0]


net <- graph_from_incidence_matrix(as.matrix(dt))
V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
V(net)$color[grep("i", names(V(net)))] <- "gold"
V(net)$color[grep("PA(16:0/18:2)", names(V(net)))] <- "black"
V(net)$shape <- c("square", "circle")[V(net)$type+1]
tkplot(net)
plot(net, vertex.frame.color = "White", vertex.label.color = "black")
