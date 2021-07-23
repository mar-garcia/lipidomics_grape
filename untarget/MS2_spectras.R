library(Spectra)
library(OrgMassSpecR)
library(Rdisop)
library(MetaboCoreUtils)
library(openxlsx)

load("maturation/data/RData/MS2_library_POS.RData")
ms2_POS <- ms2
load("tissues/data/RData/MS2_library_POS.RData")
ms2_POS <- c(ms2_POS, ms2)
load("maturation/data/RData/MS2_library_NEG.RData")
ms2_NEG <- ms2
load("tissues/data/RData/MS2_library_NEG.RData")
ms2_NEG <- c(ms2_NEG, ms2)
rm(ms2)

ions.main <- data.frame(rbind(
  c("FFA",     "[M+H]+",   "[M+CHO2]-"),
  c("CAR",     "[M+H]+",   "[M+CHO2]-"),
  c("CER",     "[M+H]+",   "[M+CHO2]-"),
  c("PA",      "[M+NH4]+", "[M-H]-"),
  c("mPA",     "[M+NH4]+", "[M-H]-"),
  c("dmPA",    "[M+NH4]+", "[M-H]-"),
  c("Lyso_PE", "[M+H]+",   "[M-H]-"),
  c("PE",      "[M+H]+",   "[M-H]-"),
  c("Lyso_PC", "[M+H]+",   "[M+CHO2]-"),
  c("PC",      "[M+H]+",   "[M+CHO2]-"),
  c("PG",      "[M+NH4]+", "[M-H]-"),
  c("PI",      "[M+NH4]+", "[M-H]-"),
  c("PS",      "[M+H]+", "[M-H]-"),
  c("MGDG",    "[M+Na]+", "[M+CHO2]-"),
  c("DGDG",    "[M+NH4]+", "[M+CHO2]-"),
  c("DAG",     "[M+NH4]+", "[M+CHO2]-"),
  c("TAG",     "[M+NH4]+", "[M+CHO2]-"),
  c("SM",      "[M+H]+",   "[M-H]-"),
  c("unknown", "[M+H]+",   "[M-H]-"),
  c("others", "[M+H]+",   "[M-H]-")
))
colnames(ions.main) <- c("class", "POS", "NEG")


data <- read.csv("output/ids_table.csv")
ms2_dt <- read.xlsx("output/MS2_spectra.xlsx")
data <- data[!data$compound %in% ms2_dt$compound,]
for(i in seq(nrow(data))){
  for(p in c("POS", "NEG")){
    if(grepl("D", data$formula[i])){
      i.mass <- getMolecule(data$formula[i])$exactmass
    } else if(!grepl("C", data$formula[i])){
      i.mass <- as.numeric(data$formula[i])
    } else {
      i.mass <- MonoisotopicMass(formula = ListFormula(data$formula[i]))
    }
    i.mz <- as.numeric(mass2mz(i.mass, ions.main[ions.main$class == data$class[i], p]))
    i.ms2 <- get(paste0("ms2_", p))
    i.ms2 <- filterPrecursorMz(i.ms2, i.mz + 0.01 * c(-1, 1))
    i.ms2 <- filterRt(i.ms2, mean(c(data$RT_tissues[i], data$RT_maturation[i]), na.rm = T)*60 + 5 * c(-1, 1))
    if(length(i.ms2) > 0){
      if(length(i.ms2) > 5){
        i.ms2 <- i.ms2[order(precursorIntensity(i.ms2), decreasing = T)[1:5]]
      }
      i.ms2_dt <- data.frame(
        compound = rep(data$compound[i], length(i.ms2)),
        polarity = rep(p, length(i.ms2)),
        adduct = rep(ions.main[ions.main$class == data$class[i], p]),
        file = basename(dataOrigin(i.ms2)),
        scan = scanIndex(i.ms2)
      )
      ms2_dt <- rbind(ms2_dt, i.ms2_dt)
    }
  }
}
ms2_dt <- ms2_dt[order(ms2_dt$compound, ms2_dt$polarity, ms2_dt$adduct, 
                       ms2_dt$file, ms2_dt$scan),]
write.xlsx(ms2_dt, "output/MS2_spectra.xlsx", row.names = FALSE)
