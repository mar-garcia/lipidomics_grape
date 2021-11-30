ms2_sim_neg <- function(exp, thr, prec){
  thrx <- thr[unlist(matchWithPpm(prec, 
                                  thr$neg_precursor, ppm = 10)),]
  if(nrow(thrx) > 0){
    out <- lapply(thrx$neg_ms2, function(x){SpectrumSimilarity(
      exp, x, print.graphic = FALSE, t = 0.1, b = 10)})
    out <- unlist(out)
    names(out) <- thrx$id
    out <- out[order(out, decreasing = TRUE)]
  } else {
    out <- NULL
  }
  out
}

precursor_neg <- function(class, formula){
  add <- case_when(
    class == "CAR" ~ "[M-H]-",
    
    class == "LPC" ~ "[M+CHO2]-",
    class == "LPC-Na" ~ "[M-H]-",
    class == "LPI" ~ "[M-H]-",
    
    class == "PA" ~ "[M-H]-",
    class == "mPA" ~ "[M-H]-",
    class == "PC" ~ "[M+CHO2]-",
    class == "PE" ~ "[M-H]-",
    class == "PG" ~ "[M-H]-",
    class == "PI" ~ "[M-H]-",
    class == "PS" ~ "[M-H]-",
    
    class == "CER" ~ "[M-H]-",
    class == "CER-FA" ~ "[M+CHO2]-",
    class == "CERhex" ~ "[M-H]-",
    class == "CERhex-FA" ~ "[M+CHO2]-",
    class == "CERlact" ~ "[M-H]-",
    
    class == "MGDG" ~ "[M+CHO2]-",
    class == "MGDG-Na" ~ "[M-H]-",
    class == "DGDG" ~ "[M+CHO2]-",
    class == "DGDG-Na" ~ "[M-H]-",
    
    class == "DAG" ~ "[M-H]-",
    class == "TAG" ~ "[M-H]-"
  ) 
  mass2mz(MonoisotopicMass(formula = ListFormula(formula)), add)
}

ms2_neg_std <- function(class, formula, sn = NULL){
  if(class == "CAR"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "LPC"){
    out <- list(cbind(
      mz = MonoisotopicMass(formula = ListFormula(subtractElements(formula, "CH3"))), 
      intensity = 100))
  } else if(class == "LPC-Na"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "LPI"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [M - H - inositol - H2O]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C6H12O6"))), "[M-H]-")
    
    # [glycerophosphoinositol - H - H2O]-
    ix <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C9H19O11P", "H2O"))), "[M-H]-")
    
    # [phosphoinositol - H - 2(H2O)]-
    iz <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("H3PO4C6H12O6", "H2OH2O"))), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, ix, iz), intensity = c(100, 50, 35, 35)))
  } else if(class == "PA"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H - CO]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "CO"))), "[M-H]-")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i3, i5, i7), intensity = c(85, 60, 100, 0.5)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6), intensity = c(75, 35, 10, 70, 30, 100)))
    }
  } else if(class == "mPA"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "PC"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    #  [M - CH3]-
    i3 <-  MonoisotopicMass(formula = ListFormula(subtractElements(formula, "CH3")))
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i2, i3), intensity = c(1, 100)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3), intensity = c(0.5, 0.5, 100)))
    }
  } else if(class == "PE"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [sn2 - H - H2O]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(sn2_fml, "H2O"))), "[M-H]-")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i3, i5, i7), intensity = c(100, 15, 2, 0.5)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6), 
                        intensity = c(50, 100, 10, 20, 2, 2)))
    }
  } else if(class == "PG"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H - sn1 - glycerol]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn1_fml), "C3H8O3"), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 - glycerol]-
    i8 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn2_fml), "C3H8O3"), "H2O"))), "[M-H]-")
    
    # [M-H-glycerol]-
    i9 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, "C3H8O3"), "H2O"))), "[M-H]-")
    
    # [sn2 - H - H2O]-
    i10 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(
      sn2_fml, "H2O"))), "[M-H]-")
    
    # [glycerophosphoglycerol-H-H2O]-
    i11 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C6H15O8P", "H2O"))), "[M-H]-")
    
    #[M-H-sn-glycerol]-
    i12 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn1_fml), "C3H8O3"), "H2OH2O"))), "[M-H]-")
    
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i2, i4, i6, i8, i9, i10, i12), 
                        intensity = c(100, 15, 10, 15, 0.5, 0.5, 0.5)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11), 
                        intensity = c(50, 100, 5, 15, 5, 10, 5, 20, 1, 0.5, 0.5)))
    }
  } else if(class == "PI"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn1 - H2O]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H - sn1 - H2O - inositol]-
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn1_fml), "C6H12O6"), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 - H2O - inositol]-
    i8 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn2_fml), "C6H12O6"), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 - inositol]-
    i9 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(subtractElements(
      subtractElements(formula, sn2_fml), "C6H12O6"), "H2OH2O"))), "[M-H]-")
    
    # [glycerophosphoinositol - H - H2O]-
    ix <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C9H19O11P", "H2O"))), "[M-H]-")
    
    # [glycerophosphoinositol - H - 2(H2O)]-
    iy <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("C9H19O11P", "H2OH2O"))), "[M-H]-")
    
    # [phosphoinositol - H - 2(H2O)]-
    iz <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements("H3PO4C6H12O6", "H2OH2O"))), "[M-H]-")
    
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i3, i5, i7, i9, iy, iz), 
                        intensity = c(75, 20, 100, 70, 5, 10, 5)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i7, i8, i9, ix, iy, iz), 
                        intensity = c(60, 20, 2, 20, 20, 100, 10, 70, 5, 2, 15, 10)))
    }
  } else if(class == "PS"){
    # [M - H - serine]- 100%
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      gsub("N", "", subtractElements(formula, "C3H5O2")))), "[M-H]-")
    
    # [M - H - serine - sn1]- 10%
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(gsub("N", "", subtractElements(formula, "C3H5O2")), 
                                   sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - serine - sn1 - H2O]- 20%
    sn1 <- sn[1]
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(gsub("N", "", subtractElements(formula, "C3H5O2")), 
                       sn1_fml))), "[M-H]-")
    
    # [sn1 - H]- 10%
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4), intensity = c(100, 10, 20, 10)))
  } else if(class == "CER"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn1 - H - H2O]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(sn1_fml, "H2O"))), "[M-H]-")
    
    # [sn2 - H - H2O]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(sn2_fml, "H2O"))), "[M-H]-")
    
    # [sn1 + N - O]-
    i4 <- MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(sn1_fml, "O"), "N")))
    
    # [sn1 - H + C2H3N]-
    i5 <-  mass2mz(MonoisotopicMass(formula = ListFormula(addElements(sn1_fml, "C2H3N"))), "[M-H]-")
    
    # [sn1 - H + C2H3N - O]-
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(addElements(sn1_fml, "C2H3N"), "O"))), "[M-H]-")
    
    # [sn2 - C2H5O]-
    i8 <- MonoisotopicMass(formula = ListFormula(subtractElements(sn2_fml, "C2H5O")))
    
    # [M - H - H2O]-
    i9 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "H2O"))), "[M-H]-")
    
    # [M - H - CH2O]-
    i10 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "CH2O"))), "[M-H]-")
    
    # [M - H - CH3OH]-
    i11 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "CH3OH"))), "[M-H]-")
    
    # [M - CH4O2]-
    i12 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "CH4O2"))), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i8, i9, i10, i11, i12), 
                      intensity = c(30, 20, 15, 10, 20, 80, 5, 10, 100, 100, 35)))
    
  } else if(class == "CER-FA"){
    out <- list(cbind(
      mz = mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-"), 
      intensity = 100))
  } else if(class == "CERhex"){
    out <- list(cbind(
      mz = mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(
        formula, "C6H10O5"))), "[M-H]-"),
      intensity = 100))
  } else if(class == "CERhex-FA"){
    out <- list(cbind(
      mz = mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-"), 
      intensity = 100))
  } else if(class == "CERlact"){
    # [M - H - hexose]-
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C6H10O5"))), "[M-H]-")
    
    # [M - H - hexose - H2O]-
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C6H12O6"))), "[M-H]-")
    
    # [M - H - 2(hexose)]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C12H20O10"))), "[M-H]-")
    
    out <- list(cbind(mz = c(i1, i2, i3), intensity = c(100, 65, 75)))
    
  } else if(class == "MGDG"){
    # [sn1 - H]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn1_fml)), "[M-H]-")
    
    # [sn2 - H]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(sn2_fml)), "[M-H]-")
    
    # [M - H - sn1]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    # [M - H ]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i3, i5), intensity = c(100, 25, 55)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5), intensity = c(25, 75, 10, 15, 100)))
    }
  } else if(class == "MGDG-Na"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "DGDG"){
    # [M - H - sn1 - H2O]-
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, sn1_fml))), "[M-H]-")
    
    # [M - H - sn2 - H2O]-
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, sn2_fml))), "[M-H]-")
    
    # [M - H ]-
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M-H]-")
    
    # [M - H - sn1]-
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn1_fml), "H2O"))), "[M-H]-")
    
    # [M - H - sn2 ]-
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(formula, sn2_fml), "H2O"))), "[M-H]-")
    
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i3, i4), intensity = c(25, 100, 5)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5), intensity = c(20, 15, 100, 5, 5)))
    }
  } else if(class == "DGDG-Na"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "DAG"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "TAG"){
    out <- list(cbind(mz = 0, intensity = 100))
  } 
  out[[1]]
}
