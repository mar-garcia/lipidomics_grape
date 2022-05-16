library(MetaboCoreUtils)

nl <- data.frame(rbind(
  c("[M+H]+", 1.007276, TRUE),
  c("[M+NH4]+", calculateMass("NH4")*1, TRUE),
  c("[M+H-H2O]+", (calculateMass("H2O")*(-1)) + 1.007276, TRUE),
  c("[M+H-C2H2]+", (calculateMass("C2H2")*(-1)) + 1.007276, TRUE),
  c("[M+H-hexose]+", (calculateMass("C6H10O5")*(-1)) + 1.007276, TRUE),
  c("[M+H-hexose-H2O]+", (calculateMass("C6H12O6")*(-1)) + 1.007276, TRUE),
  c("[M+H-2(hexose)]+", (calculateMass("C6H10O5C6H10O5")*(-1)) + 1.007276, TRUE),
  c("[M+H-2(hexose)-H2O]+", (calculateMass("C6H10O5C6H12O6")*(-1)) + 1.007276, TRUE),
  c("[M+H-2(hexose)-2(H2O)]+", (calculateMass("C6H12O6C6H12O6")*(-1)) + 1.007276, TRUE),
  c("[M+H-PA]+", (calculateMass("H3PO4")*(-1)) + 1.007276, TRUE),
  c("[M+H-mPA]+", (calculateMass(addElements("H3PO4", "CH2"))*(-1)) + 1.007276, TRUE),
  c("[M+H-dmPA]+", (calculateMass(addElements("H3PO4", "CH2CH2"))*(-1)) + 1.007276, TRUE),
  c("[M+H-PA-C2H2]+", (calculateMass(addElements("H3PO4", "C2H2"))*(-1)) + 1.007276, TRUE),
  c("[M+H-157.0504]+", 157.0504*(-1) + 1.007276, TRUE),
  c("[M+H-PE]+", (calculateMass("C2H8NO4P")*(-1)) + 1.007276, TRUE),
  c("[M+NH4-C3H9O6]+", (calculateMass("C3H9O6")*(-1)) + calculateMass("NH4"), TRUE),
  c("[M+H-PI]+", (calculateMass("C6H13O9P")*(-1)) + 1.007276, TRUE),
  c("[M+H-PS]+", (calculateMass("C3H8NO6P")*(-1)) + 1.007276, TRUE),
  c("[M+H-SQ]+", (calculateMass("C6H8O4SO3H2")*(-1)) + 1.007276, TRUE),
  c("[M+H-SQ-H2O]+", (calculateMass("C6H10O5SO3H2")*(-1)) + 1.007276, TRUE),
  
  c("[M-H]-",  - 1.007276, FALSE),
  c("[M+CHO2]-", calculateMass("CHO2")*1, FALSE),
  c("[M-CH3]-", (calculateMass("CH3")*(-1)), FALSE),
  c("[M-H-H2O]-", (calculateMass("H2O")*(-1)) - 1.007276, FALSE),
  c("[M-H-CO2]-", (calculateMass("CO2")*(-1)) - 1.007276, FALSE),
  c("[M-H-serine]-", (calculateMass("C3H5NO2")*(-1)) - 1.007276, FALSE),
  c("[M-H-hexose]-", (calculateMass("C6H10O5")*(-1)) - 1.007276, FALSE),
  c("[M-H-hexose-H2O]-", (calculateMass("C6H12O6")*(-1)) - 1.007276, FALSE),
  c("[M-H-hexose-2(H2O)]-", (calculateMass("C6H12O6H2O")*(-1)) - 1.007276, FALSE),
  c("[M-H-hexose-2(H2O)-CH2O]-", (calculateMass("C6H12O6H2OCH2O")*(-1)) - 1.007276, FALSE),
  c("[M-H-CH4O2]-", (calculateMass("CH4O2")*(-1)) - 1.007276, FALSE)
))

# add losses of fatty-acyl chains:
C <- seq(from = 10, to = 26, by = 1)
db <- seq(from = 0, to = 3, by = 1)
sn <- paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
            sep = ":")
sn <- cbind(sn, 
            paste0("C", gsub(":.*", "", sn), 
                   "H", as.numeric(gsub(":.*", "", sn))*2 - 
                     2*as.numeric(gsub(".*:", "", sn)), 
                   "O2"))
colnames(sn) <- c("sn", "formula")
for(i in seq(nrow(sn))){
  nl <- rbind(
    nl,
    # [M+H-18:2-H2O]+
    c(paste0("[M+H-", sn[i, "sn"],"-H2O]+"), 
      calculateMass(sn[i, "formula"])*(-1) + 1.007276, TRUE),
    # [M-H-18:2-H2O]-
    c(paste0("[M-H-", sn[i, "sn"],"-H2O]-"), 
      calculateMass(sn[i, "formula"])*(-1) - 1.007276, FALSE),
    # [M-H-18:2-H2O-glycerol]-
    c(paste0("[M-H-", sn[i, "sn"],"-H2O-glycerol]-"), 
      calculateMass(subtractElements(addElements(sn[i, "formula"], "C3H8O3"),
                                   "H2O"))*(-1) - 1.007276, FALSE),
    # [M-H-18:2-H2O-inositol]-
    c(paste0("[M-H-", sn[i, "sn"],"-H2O-inositol]-"), 
      calculateMass(subtractElements(addElements(sn[i, "formula"], "C6H12O6"),
                                   "H2O"))*(-1) - 1.007276, FALSE),
    # [M-H-18:2-H2O-serine]-
    c(paste0("[M-H-", sn[i, "sn"],"-H2O-serine]-"), 
      calculateMass(addElements(sn[i, "formula"],"C3H5NO2"))*(-1) - 1.007276, FALSE),
    # [M+H-18:2]+
    c(paste0("[M+H-", sn[i, "sn"],"]+"), 
      calculateMass(subtractElements(sn[i, "formula"], 
                                   "H2O"))*(-1) + 1.007276, TRUE),
    # [M-H-18:2]-
    c(paste0("[M-H-", sn[i, "sn"],"]-"), 
      calculateMass(subtractElements(sn[i, "formula"], 
                                   "H2O"))*(-1) - 1.007276, FALSE)
  )}
colnames(nl) <- c("name", "mass_add", "positive")
nl$mass_add <- as.numeric(nl$mass_add)
nl$mass_multi <- 1
nl$charge <- 1
nl$formula_add <- "X"
nl$formula_sub <- "X"
nl <- nl[, c("name", "mass_multi", "mass_add", "formula_add", "formula_sub", 
             "charge", "positive")]
rownames(nl) <- nl$name


ions <- data.frame(rbind(
  c(mass2mz(calculateMass("C5H14NO4P"), "[M+H]+"), "[PC+H]+", TRUE),
  c(mass2mz(calculateMass("C29H50O"), "[M+H-H2O]+"), "[ST+H-H2O]+", TRUE)
))
colnames(ions) <- c("mz", "name", "positive")
for(i in seq(nrow(sn))){
  ions <- rbind(
    ions,
    c(mass2mz(calculateMass(sn[i, "formula"]), "[M-H]-"), 
      paste0("[",sn[i, "sn"], "-H]-"), FALSE),
    c(mass2mz(calculateMass(addElements(sn[i, "formula"], "C3H4O")), "[M+H]+"), 
      paste0("[",sn[i, "sn"], "+H+C3H4O]+"), TRUE))
}
C <- seq(from = 32, to = 40, by = 1)
db <- seq(from = 0, to = 6, by = 1)
sn <- paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
            sep = ":")
sn <- cbind(sn, 
            paste0("C", gsub(":.*", "", sn), 
                   "H", as.numeric(gsub(":.*", "", sn))*2 - 
                     2*as.numeric(gsub(".*:", "", sn)), 
                   "O2"))
colnames(sn) <- c("sn", "formula")
for(i in seq(nrow(sn))){
  ions <- rbind(
    ions,
    c(mass2mz(calculateMass(addElements(sn[i, "formula"], "C3H2O2")), "[M+H]+"), 
      paste0("[DAG(", sn[i, "sn"], ")+H-H2O]+"), TRUE),
    c(mass2mz(calculateMass(addElements(sn[i, "formula"], "C3H4O3")), "[M+H]+"), 
      paste0("[DAG(", sn[i, "sn"], ")+H]+"), TRUE)
  )
}
ions$mz <- as.numeric(ions$mz)

ms2_ann <- function(mz, i, add, nl){
  sps <- data.frame(cbind(mz = mz, i = i))
  sps$i <- sps$i/max(sps$i)
  sps$ann <- NA
  for(j in seq(nrow(nl))){
    idx <- which(MsCoreUtils::between(
      sps$mz,
      as.numeric(mass2mz(mz2mass(i_mz, add), nl[j,])) + 0.01 * c(-1,1)
    ))
    if(length(idx) > 0){
      sps$ann[idx] <- paste(sps$ann[idx], "\n", rownames(nl)[j])
    }
  }
  
  ins <- matchWithPpm(sps$mz, tmp_ions$mz, ppm = 10)
  names(ins) <- seq(nrow(sps))
  ins <- ins[lapply(ins,length)>0]
  if(length(ins) > 0){
    for(j in seq(length(ins))){
      sps$ann[as.numeric(names(ins)[j])] <- paste(
        sps$ann[as.numeric(names(ins)[j])], "\n", tmp_ions$name[ins[[j]]])
    }
  }
  
  sps$ann <- gsub("NA \n ", "", sps$ann)
  return(sps)
}

rm(C, db, i, sn)
save.image("data/MS2_annotation.RData")