library(xcms)
raw_data <- readMSData(files = "C:/Users/lenovo/Documents/GitHub/lipidomics/untarget/tissues/data/x032_lipidgrape_tissues_pt11_entire_rep2_POS_FS.mzXML", mode = "onDisk")
dt <- data.frame(matrix(nrow = length(raw_data), ncol = 3))
colnames(dt) <- c("rt", "mz", "intensity")
for(i in seq(nrow(dt))){
  dt$rt[i] <- rtime(raw_data[i])/60
  idx <- which.max(unlist(intensity(raw_data[i])))
  dt$mz[i] <- unlist(mz(raw_data[i]))[idx]
  dt$intensity[i] <- unlist(intensity(raw_data[i]))[idx]
}
write.csv(as.matrix(dt), "maxmz.csv", row.names = F)
