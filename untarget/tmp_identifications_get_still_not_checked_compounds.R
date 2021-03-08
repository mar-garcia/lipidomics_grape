DT::datatable(
  features[features$C %in% cmps$C[cmps$RT_checked == 0],
           c("mzmed", "rtmed", "C", "compound", "annotation", "ppm")], 
  filter = "top")


#############################################################################

cmpsx <- cmps$C[cmps$RT_checked == 0]
tmpx <- data.frame(matrix(nrow = 0, ncol = ncol(features)))
colnames(tmpx) <- colnames(features)
for(i in 1:length(cmpsx)){
  if(length(grep(cmpsx[i], features$C)) > 0){
    tmpx <- rbind(tmpx, features[grep(cmpsx[i], features$C),])
  }
}
DT::datatable(
  tmpx[,c("mzmed", "rtmed", "C", "compound", "annotation", "ppm")], 
  filter = "top")
