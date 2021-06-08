library(CompoundDb)
feat$rtmed <- feat$rtmed/60
feat[unlist(matchWithPpm(701.5127, feat$mzmed, ppm = 5)),]

feat_grape$rtmed[unlist(matchWithPpm(701.5127, feat_grape$mzmed, ppm = 5))]/60

x <- c(18,18,16,16,14,12,12,10)
y <- c(1,2,0,1,0,0,1,0)
(idx <- which(colSums(combn(x, 3)) == 46 & colSums(combn(y, 3)) == 2))
combn(x, 3)[,idx]
combn(y, 3)[,idx]


md$coefficients["(Intercept)"] + (59*md$coefficients["C[-idx]"]) + (8*md$coefficients["db[-idx]"])
                                                                    


library(Rdisop)
#i.add <- "[M+CHO2]-"
i.add <- "[M-H]-"
i.add <- "[M+NH4]+"
for(j in 36){
  C <- j
  for(i in 0:20){
    i.mol <- paste0("C", C+3, "H", C*2 - (2 + 2*i) + 7,"O8P") # PA
    #i.mol <- paste0("C", C+4, "H", C*2 - (2 + 2*i) + 9,"O8P") # mPA
    #i.mol <- paste0("C", C+3, "H", C*2 - (2 + 2*i) + 6,"O5") # DAG
    i.mz <- unlist(mass2mz(getMolecule(i.mol)$exactmass, i.add))
    idx <- unlist(matchWithPpm(i.mz, feat$mzmed, ppm = 5))
    if(length(idx)>0){
      print(i)
      print(i.mol)
      print(feat[idx,])
    } 
  }
}




unlist(mass2mz(getMolecule("C39H68O5")$exactmass, adduct = 
                 c("[M+H-H2O]+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[2M+H]+", "[2M+Na]+",
                   "[M-H]-", "[M-H+HCOOH]-", "[M+Cl]-", "[2M-H]-")))
unlist(mass2mz(getMolecule("C39H68O5")$exactmass, "[2M+H]+")) + 17.02654
unlist(mass2mz(getMolecule("C39H68O5")$exactmass, "[2M-H]-")) + 46.00547

unlist(matchWithPpm(1278.0109, unlist(mass2mz(getMolecule("C39H68O5")$exactmass, adduct = adducts())), ppm = 10))

data <- read.csv("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/compounds_target.csv")
data <- data[data$class == "DAG", ]
plot(data$RT, data$C, col = "white")
text(data$RT, data$C, #paste(gsub("DAG", "", data$ID), "\n", 
     data$db)#)
idx <- which(data$formula == "C39H70O5")
text(data$RT[idx], data$C[idx], data$db[idx], col = "red")
data[idx,]

