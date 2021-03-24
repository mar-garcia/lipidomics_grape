library(CluMSID)
library(Rdisop)
library(CompoundDb)
library(xcms)
library(MsCoreUtils)

polarity <- "POS" # specify "POS" or "NEG"
load(paste0("maturation/data/RData/MS2_library_", polarity, ".RData"))
ms2list1 <- ms2list
load(paste0("tissues/data/RData/MS2_library_", polarity, ".RData"))
ms2list <- c(ms2list1, ms2list)
rm(ms2list1)

cmps <- read.csv("compounds_RT_formula.csv")
cmps <- cmps[cmps$class == "IS", ]

xdata <- readMSData(
  files = paste0("tissues/data/", polarity, 
                 "_FS_fixed/x004_lipidgrape_tissues_xx00_blank_rep2_", 
                 polarity, "_FS.mzData"),
  mode = "onDisk")

#xdata <- readMSData(
#  files = paste0("tissues/data/", polarity, 
#                 "_FS_fixed/x006_lipidgrape_tissues_xx00_STDmix_rep2_", 
#                 polarity, "_FS.mzData"),
#  mode = "onDisk")



k <- 7
cmps$ID[k]
(rt <- cmps$RT[k]*60)
if(polarity == "POS"){
  k.add <- "[M+H]+"
} else if(polarity == "NEG"){
  k.add <- "[M-H]-"
}
mz <- unlist(mass2mz(getMolecule(cmps$formula[k])$exactmass, k.add)) 
chr <- chromatogram(xdata, mz = mz + 0.1 * c(-1, 1))
par(mfrow=c(1,1))
plot(chr)
abline(v=rt, col = "red")
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
rt <- 946.568       
par(mfrow=c(1,2), mar = c(4,2,2,1))
plot(chr, xlim = c(rt - 50, rt + 50))
abline(v = rt, lty = 2, col = "grey")
sps <- as.data.frame(xdata[[closest(rt, rtime(xdata))]])
plot(sps$mz, sps$i, type = "h", #xlim = c(mz - 10, mz + 30), 
     xlab = "m/z", ylab = "intensity", 
     main = paste0("FS at ", sprintf("%.2f", round(rt/60, 2)), "' (", round(rt), "'')"))
idx <- which((sps$i/max(sps$i))*100 > 30)
if(polarity == "POS"){
  mzadd <- matchWithPpm(unlist(mass2mz(getMolecule(cmps$formula[k])$exactmass, 
                                       c("[M+H]+", "[M+NH4]+", "[M+Na]+"))), sps$mz, ppm = 10)
} else if(polarity == "NEG"){
  mzadd <- matchWithPpm(unlist(mass2mz(getMolecule(cmps$formula[k])$exactmass, 
                                       c("[M-H]-", "[M-H+HCOOH]-"))), sps$mz, ppm = 10)
}
for(i in seq(length(mzadd))){
  if(length(mzadd[[i]]) == 0){
    mzadd[[i]] <- NA
  }
  points(sps$mz[mzadd[[i]]], sps$i[mzadd[[i]]], type = "h", col = i+1)
  text(sps$mz[mzadd[[i]]], sps$i[mzadd[[i]]], round(sps$mz[mzadd[[i]]], 4), col = i+1, cex = 0.8)
}
idx <- idx[!idx %in% unlist(mzadd)]
if(length(idx) > 0){
  text(sps$mz[idx], sps$i[idx], round(sps$mz[idx], 4), cex = 0.8)
}

###################################

if(polarity == "POS"){
  unlist(mass2mz(getMolecule(cmps$formula[k])$exactmass, c("[M+H]+", "[M+NH4]+", "[M+Na]+")))
} else if(polarity == "NEG"){
  unlist(mass2mz(getMolecule(cmps$formula[k])$exactmass, c("[M-H]-", "[M-H+HCOOH]-")))
}
mz <- 828.5625              
ms2sub <- getSpectrum(ms2list, "precursor", mz, mz.tol = 0.01)
ms2sub <- getSpectrum(ms2sub, "rt", rt, rt.tol = 15)
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
      files = paste0(study, "/data/", polarity, "_DDA_mzmL/", ms2sub[[j]]@annotation), 
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
      files = paste0(study, "/data/", polarity, "_DDA_mzmL/", ms2sub[[j]]@annotation), 
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

