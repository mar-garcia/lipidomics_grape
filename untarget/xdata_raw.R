library(xcms)
for(s in c("tissues", "maturation")){
  for(p in c("POS", "NEG")){
    xdata <- list.files(paste0(s, "/data/", p, "_FS_fixed/"))
    if(s == "tissues"){
      xdata <- xdata[grep(paste(c("seeds", "pulp", "skin", "entire", "QC"), 
                                collapse = "|"), xdata)]
    } else if(s == "maturation"){
      xdata <- xdata[grep(paste(c("Pt", "QC"), collapse = "|"), xdata)]
      xdata <- xdata[!grepl(("MIX"), xdata)]
    }
    xdata <- xdata[!grepl("QCeq", xdata)]
    xdata <- readMSData(paste0(s, "/data/", p, "_FS_fixed/", xdata), 
                        mode = "onDisk")
    save(xdata, file = paste0("output/xdata_raw_", s, "_",  p, ".RData"))
  }
}