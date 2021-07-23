library(OrgMassSpecR)
library(Rdisop)
library(MetaboCoreUtils)
library(xcms)
library(plotly)
library(Spectra)
.ppm <- function(x, ppm = 10) {
  ppm * x / 1e6
}
matchWithPpm <- function(x, y, ppm = 0) {
  lapply(x, function(z, ppm) {
    which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
  }, ppm = force(ppm))
}
norm_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}

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
  c("MGDG",    "[M+NH4]+", "[M+CHO2]-"),
  c("DGDG",    "[M+NH4]+", "[M+CHO2]-"),
  c("DAG",     "[M+NH4]+", "[M+CHO2]-"),
  c("TAG",     "[M+NH4]+", "[M+CHO2]-"),
  c("SM",      "[M+H]+",   "[M-H]-"),
  c("unknown", "[M+H]+",   "[M-H]-"),
  c("others", "[M+H]+",   "[M-H]-")
))
colnames(ions.main) <- c("class", "POS", "NEG")

data <- read.csv("output/ids_table.csv")
for(s in c("tissues", "maturation")){
  for(p in c("POS", "NEG")){
    xdata <- readMSData(
      paste0(s, "/data/", p, "_FS_fixed/",
             list.files(paste0(s, "/data/", p, "_FS_fixed/"))[grep(
               "QC_rep6|xQC_10", list.files(paste0(s, "/data/", p, "_FS_fixed/")
               ))]), mode = "onDisk")
    assign(paste0("xdata_", s, "_", p), xdata)
  }
}
load("output/feat.RData")

load("maturation/data/RData/MS2_library_POS.RData")
ms2_POS <- ms2
load("tissues/data/RData/MS2_library_POS.RData")
ms2_POS <- c(ms2_POS, ms2)
load("maturation/data/RData/MS2_library_NEG.RData")
ms2_NEG <- ms2
load("tissues/data/RData/MS2_library_NEG.RData")
ms2_NEG <- c(ms2_NEG, ms2)
rm(ms2)

cmb <- combn(c("POS", "NEG", "tissues", "maturation"), 2)[,-c(1,6)]

ms2_dt <- openxlsx::read.xlsx("output/MS2_spectra.xlsx")

#i <- grep("PA34:1_FA16:0_18:1", data$compound)
i <- 89
#for(i in 6:(nrow(data))){
  rmarkdown::render(
    './id_report.Rmd',  # file 2
    output_file =  paste(gsub(":", "_", data$compound[i]), ".html", sep=''), 
    output_dir = './output/id_reports/')
#}
