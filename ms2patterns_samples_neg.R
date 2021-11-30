ms2_neg_smp <- function(class, formula, sn = NULL){
  if(class == "mPA"){
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")

    
    # [sn1 - H]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H - H2O]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(sn2_fml, "H2O"))), "[M-H]-")
    
    # [sn2 - H]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H - sn2]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn1]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i7), 
                      intensity = c(33, 1, 100, 1, 18, 1, 2)))
  } else if(class == "PE"){
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    
    # [sn1 - H]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn2]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4), 
                      intensity = c(37, 100, 24, 2)))
  } else if(class == "PI"){
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    
    # [glycerophosphoinositol - H - 2(H2O)]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C9H19O11P", "H2OH2O"))), "[M-H]-")
    
    # [glycerophosphoinositol - H - H2O]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C9H19O11P", "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn1 - inositol - H20]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(subtractElements(formula, sn1_fml), 
                                   "C6H12O6"), "H2O"))), "[M-H]-")
    
    # [M - H - sn1]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H - sn2 - inositol - H20]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(subtractElements(formula, sn2_fml), 
                                   "C6H12O6"), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 - inositol]-
    i8 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(subtractElements(formula, sn2_fml), 
                                   "C6H12O6"), "H2OH2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i9 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [phosphoinositol - H - H2O]-
    i10 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("H3PO4C6H12O6", "H2OH2O"))), "[M-H]-")
    
    # [phosphoinositol - H]-
    i11 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("H3PO4C6H12O6", "H2O"))), "[M-H]-")
    
    # [sn1 - H]-
    i12 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    i13 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13), 
                      intensity = c(15, 2, 22, 11, 3, 100, 74, 4, 21, 17, 1, 53, 25)))
    
  } else if(class == "MGDG"){
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    
    # [sn1 - H]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H ]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3), 
                      intensity = c(100, 27, 80)))
    
  }else if(class == "DGDG"){
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    
    # [M - H - sn1 - H2O]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn1]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H ]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3), 
                      intensity = c(24, 5, 100)))
    
  } else{
    out <- list(cbind(mz = 0, intensity = 100))
  }
  out[[1]]
}