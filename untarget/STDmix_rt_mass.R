library(readxl)
library(Rdisop)
library(OrgMassSpecR)

cmps <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/db_compounds.xlsx")
cmps <- cmps[which(cmps$ID_cmp == "C0612"):nrow(cmps),]
cmps$RT <- rowMeans(cmps[,c("lipidgrape_min", "lipidgrape_max")])
cmps <- subset(cmps, select = c("compound", "RT", "formula"))
cmps$mass <- NA
for(i in seq(nrow(cmps))){
  if(grepl("D", cmps$formula[i])){
    cmps$mass[i] <- getMolecule(cmps$formula[i])$exactmass
  } else{
    cmps$mass[i] <- MonoisotopicMass(formula = ListFormula(cmps$formula[i]))
  }
}


class <- read_xlsx("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx", sheet = "FORMULE")
class <- subset(class, select = c("class", "ID"))
class_IS <- read_xlsx("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx", sheet = "class_IS")
colnames(class_IS) <- c("class", "ID")
class <- rbind(class, class_IS)
colnames(class) <- c("class", "compound")
class_TAG <- data.frame(cbind(
  rep("TAG", 3),
  cmps$compound[100:102]
))
colnames(class_TAG) <- c("class", "compound")
class <- rbind(class, class_TAG)

data <- merge(cmps, class, by = "compound")
data <- data[order(data$class),]

mycols <- colorRampPalette(c("black", "gold", "yellow green", "forest green", "light sea green",
                             "steel blue", "dark blue", "dark magenta", "tomato",
                             "dark orange", "orange", "grey"))(length(unique(data$class)))
names(mycols) <- unique(data$class)

plot(data$RT, data$mass, xlab = "RT", ylab = "mass", 
     xlim = c(1, 30), ylim = c(100, 1100),
     pch = 16, col = mycols[data$class], cex = 1.2)
abline(v = seq(0, 30, 2.5), h = seq(100, 1000, 100), lty = 3, col = "grey")
legend("topright", pch = 16, legend = names(mycols), col = mycols, ncol = 1, bty = "n", cex = 0.8)
idx <- which(data$class == "TAG")
points(data$RT[idx], data$mass[idx], pch = 8)


db <- read.csv("std_db.csv")
data <- merge(data, db, by = "compound")

plot(data$RT, as.numeric(substr(data$formula, 2, 3)), xlab = "RT", ylab = "Carbons", 
     xlim = c(1, 30), #ylim = c(100, 1100),
     pch = 16, col = "white", cex = 1.2)
text(data$RT, as.numeric(substr(data$formula, 2, 3)), data$db, 
     col = mycols[data$class])
legend("topright", pch = 16, legend = names(mycols), col = mycols, ncol = 1, bty = "n", cex = 0.8)



###########################################################


library(readxl)
library(Rdisop)
library(OrgMassSpecR)
library(xcms)
library(Spectra)
library(CompoundDb)

cmps <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/db_compounds.xlsx")
cmps <- cmps[which(cmps$ID_cmp == "C0612"):nrow(cmps),]
cmps$RT <- rowMeans(cmps[,c("lipidgrape_min", "lipidgrape_max")])
cmps <- subset(cmps, select = c("compound", "RT", "formula"))
cmps$mass <- NA
for(i in seq(nrow(cmps))){
  if(grepl("D", cmps$formula[i])){
    cmps$mass[i] <- getMolecule(cmps$formula[i])$exactmass
  } else{
    cmps$mass[i] <- MonoisotopicMass(formula = ListFormula(cmps$formula[i]))
  }
}
class <- read_xlsx("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx", sheet = "FORMULE")
class <- subset(class, select = c("class", "ID"))
class_IS <- read_xlsx("C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx", sheet = "class_IS")
colnames(class_IS) <- c("class", "ID")
class <- rbind(class, class_IS)
colnames(class) <- c("class", "compound")
class_TAG <- data.frame(cbind(
  rep("TAG", 3),
  cmps$compound[100:102]
))
colnames(class_TAG) <- c("class", "compound")
class <- rbind(class, class_TAG)
cmps <- merge(cmps, class, by = "compound")


for(p in c("POS", "NEG")){
  pth <- paste0("tissues/data/", p, "_FS_fixed/") 
  fls <- list.files(pth) 
  fls <- fls[grep("STDmix", fls)]
  xdata <- readMSData(paste0(pth , fls), mode = "onDisk")
  assign(paste0("xdata_", p), xdata)
  
  pth <- paste0("tissues/data/", p, "_DDA_mzML/") 
  fls <- list.files(pth) 
  fls <- fls[grep("STDmix", fls)]
  ms2spec <- Spectra(paste0(pth, fls), backend = MsBackendDataFrame())
  ms2spec <- ms2spec[msLevel(ms2spec) > 1]
  assign(paste0("ms2_", p), ms2spec)
}
rm(p, pth, fls, xdata, ms2spec)

#xdata_POS2 <- filterFile(xdata_POS, 2)
#xdata_NEG2 <- filterFile(xdata_NEG, 2)

cmpsx <- cmps[cmps$class == "TAG", ]
#cmpsx$POS <- NA
#cmpsx$NEG <- NA
#for(i in seq(nrow(cmpsx))){
#  cmpsx$POS[i] <- unlist(mass2mz(cmpsx$mass[i], "[M+NH4]-"))
#  cmpsx$NEG[i] <- unlist(mass2mz(cmpsx$mass[i], "[M-H]-"))
#}

#for(i in seq(nrow(cmpsx))){
#  chr <- chromatogram(xdata_POS2, mz = cmpsx$POS[i] + 0.01 * c(-1, 1))
#  assign(paste0("chr_POS_", formatC(i, width = 2, flag = "0")), chr)
#  
#  chr <- chromatogram(xdata_NEG2, mz = cmpsx$NEG[i] + 0.01 * c(-1, 1))
#  assign(paste0("chr_NEG_", formatC(i, width = 2, flag = "0")), chr)
#}

#plot(0,0, xlim = c(0, 30*60), ylim = c(0, 1.13), col = "white", xlab = "RT", ylab = "relative intensity", xaxt="n")
#axis(1, at = seq(0, 60*30, 60), labels = seq(0, 30))

#for(i in seq(nrow(cmpsx))){
#  chr <- get(paste0("chr_NEG_", formatC(i, width = 2, flag = "0")))
#  points(rtime(chr[[1]]), intensity(chr[[1]])/max(intensity(chr[[1]]), na.rm = T), type = "l", col = i)
#  text(cmpsx$RT[i]*60, 1.07, cmpsx$compound[i], srt = 90, cex = 0.5, col = i)
#}

cmpsx$RT[cmpsx$compound == "17:1 Lyso PI"] <- 5.84+0.33
cmpsx$compound[cmpsx$compound == "15:0-18:1(d7)-15:0 TAG_Na"] <- "TAG48:1_FA18:1_d7"

plot(cmpsx$RT, as.numeric(substr(cmpsx$formula, 2, 3)), 
     col = "white", xlab = "RT", ylab = "nº C", xlim = c(14.3, 22), ylim = c(26,58))
text(cmpsx$RT, as.numeric(substr(cmpsx$formula, 2, 3)),
     paste(substr(cmpsx$compound, 4, 14), "\n", 
           (substr(cmpsx$compound, 7, 7))
           ))
abline(h = seq(16,18, 2), lty = 3, col = "grey")

C <- as.numeric(substr(cmpsx$formula, 2, 3))
db <- as.numeric(substr(cmpsx$compound, 7, 7))
rt <- cmpsx$RT

md <- lm(rt~C+db)
summary(md)


rt_pred <- md$coefficients["(Intercept)"] + C*md$coefficients["C"] + db*md$coefficients["db"]
plot(rt, rt_pred)

xC <- 61
xdb <- 7
xrt <- 21.94
(md$coefficients["(Intercept)"] + xC*md$coefficients["C"] + xdb*md$coefficients["db"]) - abs(min(md$residuals))
(md$coefficients["(Intercept)"] + xC*md$coefficients["C"] + xdb*md$coefficients["db"]) + abs(max(md$residuals))



((md$coefficients["(Intercept)"] + (xC*md$coefficients["C"]) - xrt) / -md$coefficients["db"])

