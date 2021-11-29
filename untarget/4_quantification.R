library(xcms)
library(MetaboCoreUtils)
library(Rdisop)
.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}

cmps <- read.csv("compounds.csv")
res <- read.csv("output/res.csv")
res$id <- res$compound
res$id <- gsub(")", "", res$id)
res$id <- gsub("\\(", "_", res$id)
res$id <- gsub("/", "_", res$id)

for(s in c("tissues", "maturation")){
  if(!dir.exists(paste0("output/chr_", s))){
    dir.create(paste0("output/chr_", s))
  }
  for(p in c("POS", "NEG")){
    load(paste0(s, "/data/RData/data_", p, ".RData"))
    cmp <- unique(db_exp$compound[!grepl("13C", db_exp$adduct)])
    db_exp <- db_exp[db_exp$compound %in% cmp, ]
    db_exp$polarity <- p
    db_exp$study <- s
    assign(paste0("dbexp_", p), db_exp)
    feat <- as.data.frame(featureDefinitions(xdata))
    assign(paste0("feat_", s, "_", p), feat)
    ydata <- list.files(paste0(s, "/data/", p, "_FS_fixed/"))
    if(s == "tissues"){
      ydata <- ydata[grep(paste(c("seeds", "pulp", "skin", "entire", "QC"), 
                                collapse = "|"), ydata)]
    } else if(s == "maturation"){
      ydata <- ydata[grep(paste(c("Pt", "QC"), collapse = "|"), ydata)]
      ydata <- ydata[!grepl(("MIX"), ydata)]
    }
    ydata <- ydata[!grepl("QCeq", ydata)]
    ydata <- readMSData(paste0(s, "/data/", p, "_FS_fixed/", ydata), 
                        mode = "onDisk")
    assign(paste0("ydata_", s, "_", p), ydata)
    
    cwp <- processParam(processHistory(xdata)[[1]])
    int <- cwp@prefilter[2]
    intx <- as.numeric(paste0(substr(
      as.character(int), 1, nchar(as.character(int))-1), 
      as.numeric(substr(as.character(int), nchar(as.character(int)), 
                        nchar(as.character(int))))-1))
    cwp@prefilter[2] <- intx
    xchr <- findChromPeaks(ydata, param = cwp)
    cp <- chromPeaks(xchr)
    assign(paste0("cp_", s, "_", p), cp)
    assign(paste0("xdata_", s, "_", p), xdata)
    
    col_smpl <- rep("grey", length(fileNames(ydata)))
    col_smpl[grep("skin", fileNames(ydata))] <- "#377EB8"
    col_smpl[grep("pulp", fileNames(ydata))] <- "#4DAF4A"
    col_smpl[grep("seeds", fileNames(ydata))] <- "#984EA3"
    col_smpl[grep("entire", fileNames(ydata))] <- "#A65628"
    col_smpl[grep("QCdl_02uL", fileNames(ydata))] <- "#FFFFB2"
    col_smpl[grep("QCdl_05uL", fileNames(ydata))] <- "#FED976"
    col_smpl[grep("QCdl_10uL", fileNames(ydata))] <- "#FEB24C"
    col_smpl[grep("QCrw_02uL", fileNames(ydata))] <- "#FD8D3C"
    col_smpl[grep("QCrw_05uL", fileNames(ydata))] <- "#F03B20"
    col_smpl[grep("QCrw_10uL", fileNames(ydata))] <- "#BD0026"
    col_smpl[grep("Pt_01", fileNames(ydata))] <- "#EFF3FF"
    col_smpl[grep("Pt_02", fileNames(ydata))] <- "#C6DBEF"
    col_smpl[grep("Pt_03", fileNames(ydata))] <- "#9ECAE1"
    col_smpl[grep("Pt_04", fileNames(ydata))] <- "#6BAED6"
    col_smpl[grep("Pt_05", fileNames(ydata))] <- "#4292C6"
    col_smpl[grep("Pt_06", fileNames(ydata))] <- "#2171B5"
    col_smpl[grep("Pt_07", fileNames(ydata))] <- "#084594"
    col_smpl[grep("Pt_08", fileNames(ydata))] <- "#9E9AC8"
    col_smpl[grep("Pt_09", fileNames(ydata))] <- "#807DBA"
    col_smpl[grep("Pt_10", fileNames(ydata))] <- "#6A51A3"
    col_smpl[grep("Pt_11", fileNames(ydata))] <- "#54278F"
    col_smpl[grep("Pt_12", fileNames(ydata))] <- "#3F007D"
    col_smpl[grep("Pt_13", fileNames(ydata))] <- "#3F007D"
    assign(paste0("col_", s, "_", p), col_smpl)
  }
  
  dbexp <- rbind(dbexp_NEG, dbexp_POS)
  assign(paste0("dbexp_", s), dbexp)
  
  data <- matrix(nrow = length(fileNames(ydata)), ncol = nrow(res))
  if(s == "tissues"){
    smpl <- unlist(substring(basename(fileNames(ydata)), 25, 40))
    smpl <- gsub("xx00_", "", smpl)
    smpl <- gsub("_N|_P|_NEG|_POS", "", smpl)
    idx <- grepl("_$", smpl)
    smpl[idx] <- unlist(substring(smpl[idx], 1, 15))
  } else if(s == "maturation"){
    smpl <- unlist(substring(basename(fileNames(ydata)), 17, 28))
    smpl <- gsub("xx_00_xQC_", "QC", smpl)
    smpl <- gsub("Rep_", "", smpl)
    smpl <- gsub("Pt_", "Pt", smpl)
  }
  rownames(data) <- smpl
  colnames(data) <- res$id
  assign(paste0("data_", s), data)
}
dbexp <- rbind(dbexp_tissues, dbexp_maturation)
dbexp <- dbexp[!grepl("13C", dbexp$adduct),]
rm(feat, ydata, ms2, xdata, smpl, col_smpl, data, 
   db_exp, dbexp_NEG, dbexp_POS, dbexp_tissues, dbexp_maturation)

for(i in seq(nrow(res))){
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
  for(s in c("tissues", "maturation")){
    ft <- dbexp$FT[dbexp$compound == res$id[i] & dbexp$study == s]
    xdata <- get(paste0("xdata_", s, "_", pol))
    mz <- mass2mz(getMolecule(fml)$exactmass, ad)
    ydata <- get(paste0("ydata_", s, "_", pol))
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
    
    
    if(!is.na(cmps[cmps$compound == res$compound[i], paste0(substr(s, 1, 1), "_min")])){
      i.cmp[,"rtmin"] <- cmps[cmps$compound == res$compound[i], paste0(substr(s, 1, 1), "_min")]*60
    }
    if(!is.na(cmps[cmps$compound == res$compound[i], paste0(substr(s, 1, 1), "_max")])){
      i.cmp[,"rtmax"] <- cmps[cmps$compound == res$compound[i], paste0(substr(s, 1, 1), "_max")]*60
    }
    
    pks <- manualChromPeaks(ydata,
                            chromPeaks = i.cmp,
                            samples = seq_along(fileNames(ydata)),
                            BPPARAM = bpparam(),
                            msLevel = 1L)
    fl <- paste0("output/chr_", s, "/", gsub(":", "_", res$id[i]), ".jpeg")
    if(!(gsub(".*/", "", fl) %in% list.files(paste0("output/chr_", s, "/")))){
      col <- get(paste0("col_", s, "_", pol))
      jpeg(fl)
      chr <- chromatogram(ydata, 
                          mz = i.cmp[,1:2], rt = i.cmp[,3:4] + 10 * c(-1,1))
      plot(chr, col = col, xaxt = "n")
      axis(1, at = seq(0, 60*30, 3), labels = sprintf("%.2f", seq(0, 30, 3/60))) 
      abline(v=c(rtmin, rtmax), lty = 2, col = "grey")
      abline(v=i.cmp[,3:4])
      dev.off()
    }
    pks <- data.frame(chromPeaks(pks))
    if(s == "tissues"){
      tmp <- unlist(substring(basename(fileNames(ydata)), 25, 40))
      tmp <- gsub("xx00_", "", tmp)
      tmp <- gsub("_N|_P|_NEG|_POS", "", tmp)
      idx <- grepl("_$", tmp)
      tmp[idx] <- unlist(substring(tmp[idx], 1, 15))
    } else if(s == "maturation"){
      tmp <- unlist(substring(basename(fileNames(ydata)), 17, 28))
      tmp <- gsub("xx_00_xQC_", "QC", tmp)
      tmp <- gsub("Rep_", "", tmp)
      tmp <- gsub("Pt_", "Pt", tmp)
    }
    pks$smpl <- tmp[pks$sample]
    data <- get(paste0("data_", s))
    data[pks$smpl, res$id[i]] <- pks$into
    assign(paste0("data_", s), data)
  }
  print(i)
}




for(s in c("tissues", "maturation")){
  data <- get(paste0("data_", s))
  write.csv(data, paste0("output/data_", s, ".csv"))
  idx <- grep("QC", rownames(data))
  idx <- idx[!grepl("rw|dl", rownames(data)[idx])]
  CV <- (apply(data[idx,], 2, sd) / apply(data[idx,], 2, mean))*100
  assign(paste0("CV_", s), CV)
  
  idx <- grep("QC", rownames(data))
  idx <- idx[grepl("rw|dl", rownames(data)[idx])]
  dil <- rownames(data)[idx]
  dil <- gsub("uL.*", "", dil)
  dil1 <- gsub("_.*", "", dil)
  dil1 <- gsub("QCrw", 1, dil1)
  dil1 <- gsub("QCdl", 10, dil1)
  dil1 <- as.numeric(dil1)
  dil2 <- gsub(".*_", "", dil)
  dil2 <- gsub("uL", "", dil2)
  dil2 <- gsub("02", "2.5", dil2)
  dil2 <- as.numeric(dil2)
  dil <- dil2 / dil1
  rm(dil1, dil2)
  tmp <- cbind(rownames(data)[idx], dil)
  linearity <- cor(log10(data[idx, ]), log10(dil))
  assign(paste0("lin_", s), linearity)
}

res <- as.data.frame(cbind(res$compound, CV_tissues, lin_tissues, CV_maturation, lin_maturation))
colnames(res) <- c("compound", "t_CV", "t_linearity", "m_CV", "m_linearity")
res$t_CV <- sprintf("%.1f", as.numeric(res$t_CV))
res$m_CV <- sprintf("%.1f", as.numeric(res$m_CV))
res$t_linearity <- sprintf("%.3f", as.numeric(res$t_linearity))
res$m_linearity <- sprintf("%.3f", as.numeric(res$m_linearity))


DT::datatable(res, rownames = FALSE, options = list(pageLength = nrow(res)))
