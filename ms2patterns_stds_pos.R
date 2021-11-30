ms2_sim_pos <- function(exp, thr, prec){
  thrx <- thr[unlist(matchWithPpm(prec, 
                                  thr$pos_precursor, ppm = 10)),]
  if(nrow(thrx) > 0){
    out <- lapply(thrx$pos_ms2, function(x){SpectrumSimilarity(
      exp, x, print.graphic = FALSE, t = 0.1, b = 10)})
    out <- unlist(out)
    names(out) <- thrx$id
    out <- out[order(out, decreasing = TRUE)]
  } else {
    out <- NULL
  }
  out
}

precursor_pos <- function(class, formula){
  add <- case_when(
    class == "CAR" ~ "[M+H]+",
    
    class == "LPC" ~ "[M+H]+",
    class == "LPC-Na" ~ "[M+Na]+",
    class == "LPI" ~ "[M+H]+",
    
    class == "PA" ~ "[M+NH4]+",
    class == "mPA" ~ "[M+NH4]+",
    class == "PC" ~ "[M+Na]+",
    class == "PE" ~ "[M+H]+",
    class == "PG" ~ "[M+NH4]+",
    class == "PI" ~ "[M+NH4]+",
    class == "PS" ~ "[M+H]+",
    
    class == "CER" ~ "[M+H]+",
    class == "CER-FA" ~ "[M+H]+",
    class == "CERhex" ~ "[M+H]+",
    class == "CERhex-FA" ~ "[M+H]+",
    class == "CERlact" ~ "[M+H]+",
    
    class == "MGDG" ~ "[M+NH4]+",
    class == "MGDG-Na" ~ "[M+Na]+",
    class == "DGDG" ~ "[M+NH4]+",
    class == "DGDG-Na" ~ "[M+Na]+",
    
    class == "DAG" ~ "[M+NH4]+",
    class == "TAG" ~ "[M+NH4]+"
  ) 
  mass2mz(MonoisotopicMass(formula = ListFormula(formula)), add)
}


ms2_pos_std <- function(class, formula, sn = NULL){
  if(class == "CAR"){
    # [M + H - trimethylamine]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C3H9N"))), "[M+H]+")
    
    # [sn + H - H2O]+
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(sn1_fml, "H2O"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(100, 15)))
  } else if(class == "LPC"){
    # [M + H - H2O]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "H2O"))), "[M+H]+")
    
    # [Phosphocholine]+
    i2 <- MonoisotopicMass(formula = ListFormula("C5H15NO4P"))
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(100, 25)))
  } else if(class == "LPC-Na"){
    out <- list(cbind(mz = mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C3H9N"))), "[M+Na]+"), 
      intensity = 100))
  } else if(class == "LPI"){
    # [M + H - H2O]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(
      formula, "H2O"))), "[M+H]+")
    
    # [M + H - 2(H2O)]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(
      formula, "H4O2"))), "[M+H]+")
    
    # [M + H - phosphoinositol]+
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(subtractElements(formula, "C6H12O6"), "H3PO4"), 
      "H2O"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2, i3), intensity = c(100, 30, 55)))
  } else if(class == "PA"){
    # [M + H]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M+H]+")
    
    # [M+H-phosphate]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "H3PO4"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(50, 100)))
  } else if(class == "mPA"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "PC"){
    # [M + Na - trimethylamine]+ 100%
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C3H9N"))), "[M+Na]+")
    
    # [M + Na - phosphoconline]+ 95%
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C5H14NO4P"))), "[M+Na]+")
    
    # [M + H - phosphoconline]+ 5%
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C5H14NO4P"))), "[M+Na]+")
    
    # [M + Na - sn1 - H2O]+ 0.5%
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+Na]+")
    
    # [M + Na - sn2 - H2O]+ 0.1%
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M+Na]+")
    
    # [M + H - sn2 - H2O]+ 0.5%
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M+H]+")
    
    # [M + Na - sn1 - H2O - trimethylamine]+ 1%
    i7 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(subtractElements(formula, sn1_fml), "C3H9N"))), "[M+Na]+")
    
    # [M + Na - sn2 - H2O - trimethylamine]+ 0.2%
    i8 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(subtractElements(formula, sn2_fml), "C3H9N"))), "[M+Na]+")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1, i2, i3, i4, i7), 
                        intensity = c(100, 95, 5, 0.5, 12)))
    } else {
      out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6, i7, i8), 
                        intensity = c(100, 95, 5, 0.5, 0.1, 0.5, 1, 0.2)))
    }
    
  } else if(class == "PE"){
    out <- list(cbind(mz = mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C2H8NO4P"))), "[M+H]+"), 
      intensity = 100))
  } else if(class == "PG"){
    out <- list(cbind(mz = mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C3H9O6"))), "[M+NH4]+"), 
      intensity = 100))
  } else if(class == "PI"){
    # [M + H]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M+H]+")
    
    # [M + H - phosphoinositol]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(addElements(
      subtractElements(subtractElements(formula, "C6H12O6"), "H3PO4"), 
      "H2O"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(55, 100)))
  } else if(class == "PS"){
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C3H8NO6P"))), "[M+H]+")
    
    # [M + H - sn1]+ 0.3%
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(formula, sn1_fml), "H2O"))), "[M+H]+")
    
    # [M + H - sn1 - H2O]+ 0.3%
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+H]+")
    
    # [M + H - sn1 - 2(H2O)]+ 0.1%
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(subtractElements(formula, sn1_fml), "H2O"))), "[M+H]+")
    
    # [M + H - sn1 - phosphoserine]+ 0.5%
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(addElements(subtractElements(formula, sn1_fml), "H2O"), "C3H8NO6P"))), "[M+H]+")
    
    # [M + H - sn1 - H2O - phosphoserine]+ 0.2%
    i6 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(subtractElements(formula, sn1_fml), "C3H8NO6P"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5, i6), intensity = c(100, 0.3, 0.3, 0.1, 0.5, 0.2)))
    
  } else if(class == "CER"){
    # [M + H - H2O]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "H2O"))), "[M+H]+")
    
    # [M+H-2(H2O)]+ 1%
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "H4O2"))), "[M+H]+")
    
    # [sn1 + NH4]+ 1%
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      sn2_fml)), "[M+NH4]+")
    
    # [sn2 + NH4 - H2O]+ 4%
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(sn2_fml, "H2O"))), "[M+NH4]+")
    
    # [sn1 + NH4 - H2O]+ 1%
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(sn1_fml, "H2O"))), "[M+NH4]+")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5), intensity = c(100, 1, 1, 4, 1)))
    
  } else if(class == "CER-FA"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "CERhex"){
    # [M + H - H2O]+ 100%
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "H2O"))), "[M+H]+")
    
    # [M+H-hexose]+ 3%
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C6H10O5"))), "[M+H]+")
    
    # [M+H-hexose-H2O]+ 20%
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, "C6H12O6"))), "[M+H]+")
    
    # [M+H-hexose-2(H2O)]+ 2%
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(subtractElements(formula, "C6H12O6"), "H2O"))), "[M+H]+")
    
    # [sn2 + H + -HO2 + N]+ 0.5%
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i5 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      addElements(subtractElements(sn2_fml, "HO2"), "N"))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2, i3, i4, i5), intensity = c(100, 3, 20, 2, 0.5)))
  } else if(class == "CERhex-FA"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "CERlact"){
    out <- list(cbind(mz = 0, intensity = 100))
  } else if(class == "MGDG"){
    # [M + H - galactose]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C6H12O6"))), "[M+H]+")
    
    # [M + H - galactose + 0.984]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C6H12O6"))) + 0.984, "[M+NH4]+")
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(25, 100)))
  } else if(class == "MGDG-Na"){
    # [M + Na - sn1 - H2O]+
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+Na]+")
    
    # [M + Na - sn2 - H2O]+
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M+Na]+")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i1), intensity = c(100)))
    } else {
      out <- list(cbind(mz = c(i1, i2), intensity = c(45, 100)))
    }
  } else if(class == "DGDG"){
    # [M + H - 2(galactose)]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C12H22O11"))), "[M+H]+")
    
    # [M + H - 2(galactose) + 0.984]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C12H22O11"))) + 0.984, "[M+NH4]+")
    
    out <- list(cbind(mz = c(i1, i2), intensity = c(20, 100)))
  } else if(class == "DGDG-Na"){
    # [M + Na - sn1 - H2O]+
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+Na]+")
    
    # [M + Na - sn2 - H2O]+
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M+Na]+")
    
    # [M + Na - hexose]+
    #i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(formula, "C16H10O5"))), "[M+Na]+")
    
    if(sn1 == sn2){
      out <- list(cbind(mz = c(i2), intensity = c(100)))
    } else {
      out <- list(cbind(mz = c(i1, i2), intensity = c(50, 100)))
    }
  } else if(class == "DAG"){
    # [M + H]+
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M+H]+")
    
    # [M + H - H2O]+
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(subtractElements(
      formula, "H2O"))), "[M+H]+")
    
    # [M + H - sn1 - H2O]+
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+H]+")
    
    out <- list(cbind(mz = c(i1, i2, i3), intensity = c(30, 100, 65)))
  } else if(class == "TAG"){
    # [M + H]+
    i4 <- mass2mz(MonoisotopicMass(formula = ListFormula(formula)), "[M+H]+")
    
    # [M + H - sn1 - H2O]+
    sn1 <- sn[1]
    sn1_C <- as.numeric(gsub(":.*", "", sn1))
    sn1_db <- as.numeric(gsub(".*:", "", sn1))
    sn1_fml <- paste0("C", sn1_C, "H", sn1_C*2 - (2*sn1_db), "O2")
    i1 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn1_fml))), "[M+H]+")
    
    # [M + H - sn2 - H2O]+
    sn2 <- sn[2]
    sn2_C <- as.numeric(gsub(":.*", "", sn2))
    sn2_db <- as.numeric(gsub(".*:", "", sn2))
    sn2_fml <- paste0("C", sn2_C, "H", sn2_C*2 - (2*sn2_db), "O2")
    i2 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn2_fml))), "[M+H]+")
    
    # [M + H - sn3 - H2O]+
    sn3 <- sn[3]
    sn3_C <- as.numeric(gsub(":.*", "", sn3))
    sn3_db <- as.numeric(gsub(".*:", "", sn3))
    sn3_fml <- paste0("C", sn3_C, "H", sn3_C*2 - (2*sn3_db), "O2")
    i3 <- mass2mz(MonoisotopicMass(formula = ListFormula(
      subtractElements(formula, sn3_fml))), "[M+H]+")
    
    if(sn1 == sn2 & sn2 == sn3){
      out <- list(cbind(mz = c(i1), intensity = c(100)))
    } else if((sn1 != sn2 & sn2 == sn3) | sn1 != sn2 & sn1 == sn3){
      out <- list(cbind(mz = c(i1, i2), intensity = c(100, 100)))
    } else if(sn1 == sn2 & sn2 != sn3){
      out <- list(cbind(mz = c(i1, i3), intensity = c(100, 40)))
    }else {
      out <- list(cbind(mz = c(i1, i2, i3), intensity = c(100, 100, 100)))
    }
  }
  out[[1]]
}
