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

library(Spectra)
library(igraph)
library(OrgMassSpecR)
library(MetaboCoreUtils)
library(Rdisop)

# Get MS2 ----
pol <- "POS"
fls <- "C:/Users/lenovo/Documents/GitHub/lipidomics/untarget/tissues/data/POS_DDA_mzML/x009_lipidgrape_tissues_pt11_skin_rep2_POS_DDA.mzML"
ms2 <- Spectra(fls, backend = MsBackendDataFrame())
ms2 <- ms2[msLevel(ms2) > 1]
ms2 <- ms2[precursorIntensity(ms2) > 1e6]

# IDs ----
## ions ----
sn <- data.frame(
  "C" = rep(seq(14, 22), each = 4),
  "db" = rep(seq(0, 3), 9)
)
sn$formula <- paste0("C", sn$C, "H", sn$C*2 - 2*sn$db, "O2")
sn$sn <- paste0(sn$C, ":", sn$db)
sn$mass <- NA
sn.list <- list()
for(i in seq(nrow(sn))){
  sn$mass[i] <- MonoisotopicMass(formula = ListFormula(sn$formula[i]))
  
  sn.list[[i]] <- sn$sn[i]
  names(sn.list)[[i]] <- sn$sn[i]
}

ids <- readxl::read_xlsx("identification.xlsx")
ions <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(ions) <- c("class", "compound", "RT", "mz", "adduct", "polarity")
for(i in seq(nrow(ids))){
  i.ions <- unlist(strsplit(ids$ions[i], "; "))
  if(grepl("D", ids$formula[i])){
    j.mass <- Rdisop::getMolecule(ids$formula[i])$exactmass
  } else {
    j.mass <- MonoisotopicMass(formula = ListFormula(ids$formula[i]))
  }
  for(j in seq(length(i.ions))){
    j.add <- gsub(".*\\[", "[", i.ions[j])
    j.add <- gsub("-H\\+HCOOH", "\\+CHO2", j.add)
    if(j.add %in% c(adductNames(polarity = "positive"), 
                    adductNames(polarity = "negative"))){
      j.mz <- unlist(mass2mz(j.mass, j.add))
    } else if(grepl("C2H8N", j.add)){
      if(grepl("2M", j.add)){
        j.mz <- j.mass*2 + getMolecule("C2H8N")$exactmass
      } else {
        j.mz <- j.mass + getMolecule("C2H8N")$exactmass
      }
    } else if(grepl("-PA", j.add)){
      j.mz <- j.mass - getMolecule("H2PO4")$exactmass
    } else if(grepl("-PI", j.add)){
      if(grepl("NH4", j.add)){
        j.mz <- as.numeric(mass2mz(j.mass, "[M+NH4]+")) - getMolecule("H2PO4C6H12O6")$exactmass
      } else {
        j.mz <- j.mass - getMolecule("H2PO4C6H12O6")$exactmass
      }
    } else if(grepl("-mPA", j.add)){
      j.mz <- j.mass - getMolecule("H2PO4CH2")$exactmass
    } else if(grepl("-CH3", j.add)){
      j.mz <- j.mass - getMolecule("CH3")$exactmass
    } else if(grepl("-CH2", j.add)){
      if(grepl("2M", j.add)){
        j.mz <- j.mass*2 - getMolecule("CH3")$exactmass
      } else {
        j.mz <- j.mass - getMolecule("CH3")$exactmass
      }
    } else if(grepl(":", j.add)){
      if(grepl("H-", j.add)){
        j.mz <- mass2mz(j.mass,
                        paste0(substr(j.add, 1, 4), 
                               substr(j.add, nchar(j.add)-1, nchar(j.add)))) - 
          sn$mass[sn$sn == substr(j.add, 6, nchar(j.add)-2)]
      } else{
        j.mz <- sn$mass[sn$sn == gsub("\\[", "", gsub("-H]-", "", j.add))]
        j.mz <- as.numeric(
          mass2mz(j.mz, gsub(gsub("\\[", "", gsub("-H]-", "", j.add)), "M", j.add)))
      }
    }
    if(grepl("13C", i.ions[j])){
      if(grepl("\\(", i.ions[j])){
        j.mz <- j.mz + 1.003355*as.numeric(gsub("\\(.*", "", i.ions[j]))
      } else {
        j.mz <- j.mz + 1.003355
      }
    } 
    if(substr(i.ions[j], nchar(i.ions[j]), nchar(i.ions[j])) == "+"){
      i.pol <- "POS"
    } else if(substr(i.ions[j], nchar(i.ions[j]), nchar(i.ions[j])) == "-"){
      i.pol <- "NEG"
    }
    
    ions[nrow(ions)+1, ] <- c(ids$class[i], ids$compound[i], ids$RT[i], j.mz, 
                              i.ions[j], i.pol)
  }
}
ions$RT <- as.numeric(ions$RT)
ions$mz <- as.numeric(ions$mz)

## annotation ----
ms2$name <- paste(sprintf("%.4f", precursorMz(ms2)), sprintf("%.2f", rtime(ms2)/60), sep = "_")
ms2$adduct <- "x"
ions <- ions[ions$polarity == pol, ]
for(i in seq(length(ms2))){
  idx <- unlist(matchWithPpm(precursorMz(ms2[i]), ions$mz, ppm = 10))
  idx <- idx[abs(rtime(ms2[i]) - ions$RT[idx]*60) < 10]
  if(length(idx) > 0){
    ms2$name[i] <- paste(ions$compound[idx], collapse = "; ")
    ms2$adduct[i] <- paste(ions$adduct[idx], collapse =  "; ")
  } 
}

## duplicates ----
tmp <- unique(paste(ms2$name, ms2$adduct)[which(duplicated(paste(ms2$name, ms2$adduct)))])
for(i in seq(length(tmp))){
  idx <- which(ms2$name == gsub(" .*", "", tmp[i]) & ms2$adduct == gsub(".* ", "", tmp[i]))
  i.ms2 <- ms2[idx]
  i.ms2 <- addProcessing(i.ms2, norm_int)
  i.ms2 <- replaceIntensitiesBelow(i.ms2, threshold = 10, value = 0)
  i.ms2 <- filterIntensity(i.ms2, intensity = c(0.1, 100))
  if(min(compareSpectra(i.ms2) > 0.5)){
    i.ms2 <- combineSpectra(i.ms2)
    ms2 <- ms2[-idx]
    ms2 <- c(ms2, i.ms2)
  } else{
    if(length(idx) > 2){
      tmp2 <- seq(length(idx))
      for(j in seq(length(idx))){
        if(j %in% tmp2){
          j.idx <- which(compareSpectra(i.ms2[j], i.ms2) > 0.5)
          j.ms2 <- combineSpectra(i.ms2[j.idx])
          ms2 <- ms2[-idx[j.idx]]
          j.ms2$name <- paste0(j.ms2$name, letter[j])
          ms2 <- c(ms2, j.ms2)
          tmp2 <- tmp2[!tmp2 %in% j.idx]
        }
      }
    }
  }
}

ms2 <- addProcessing(ms2, norm_int)
ms2 <- replaceIntensitiesBelow(ms2, threshold = 10, value = 0)
ms2 <- filterIntensity(ms2, intensity = c(0.1, 100))


# MS2NET ----
mzdif.pos <- data.frame(rbind(
  c(MonoisotopicMass(formula = ListFormula("NH3")), "loss NH3 -> PA / mPA / PI / DAG / TAG / MGDG / DGDG"),
  c(MonoisotopicMass(formula = ListFormula("H2O")), "loss H2O -> Lyso PC"),
  c(MonoisotopicMass(formula = ListFormula("NH3H2O")), "loss NH3 & H2O -> DAG"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4NH3")), "loss NH3 & phosphate -> PA"),
  c(MonoisotopicMass(formula = ListFormula("NH3H3PO4CH2")), "loss NH3 & mPA -> mPA"),
  c(MonoisotopicMass(formula = ListFormula("C3H9N")), "loss C3H9N -> PC / Carnitine"),
  c(MonoisotopicMass(formula = ListFormula("C5H14NO4P")), "loss C5H14NO4P -> PC"),
  c(MonoisotopicMass(formula = ListFormula("H3PO4C6H12O6")) - 0.984, 
    "loss phosphoinositol (NH4 ion) -> PI"),
  c(MonoisotopicMass(formula = ListFormula("C2H8NO4P")), 
    "loss phosphoethanolamine -> PE"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> MGDG [M+Na]+")),
  cbind(sn$mass + MonoisotopicMass(formula = ListFormula("NH3")), 
        paste("loss NH3 &", sn$sn, "-> DAG / TAG")),
  
  c(MonoisotopicMass(formula = ListFormula("C6H12O6")) - 0.984, 
    "loss 1 galactose (NH4 ion) -> MGDG"),
  c(MonoisotopicMass(formula = ListFormula("C6H12O6")) + 
      MonoisotopicMass(formula = ListFormula("NH3")), 
    "loss NH4 & 1 galactose -> MGDG"),
  c(MonoisotopicMass(formula = ListFormula("C12H22O11")) - 0.984, 
    "loss 2 galactoses (NH4 ion) -> DGDG"),
  c(MonoisotopicMass(formula = ListFormula("C12H22O11")) + 
      MonoisotopicMass(formula = ListFormula("NH3")), 
    "loss NH4 & 2 galactoses -> DGDG"),
  cbind(MonoisotopicMass(formula = ListFormula("C5H15NO4P")), "[Phosphocholine]+ -> Lyso PC")
))
colnames(mzdif.pos) <- c("dif", "add")
mzdif.pos$dif <- as.numeric(mzdif.pos$dif)

mzdif.neg <- data.frame(rbind(
  c(MonoisotopicMass(formula = ListFormula("HCOOH")), "loss HCOOH -> MGDG / DGDG"),
  cbind(sn$mass, paste("loss ", sn$sn, "-> PA / PG / PI")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O")), 
        paste0("loss '", sn$sn, "-H2O' -> PA / mPA / PG / PE")),
  cbind((sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))) + 
          MonoisotopicMass(formula = ListFormula("C3H8O3")), 
        paste0("loss '", sn$sn, "-H2O' & glycerol -> PA / PG / PE")),
  cbind(sn$mass + MonoisotopicMass(formula = ListFormula("HCOOH")), 
        paste("loss HCOOH &", sn$sn, "-> DGDG")),
  cbind(sn$mass - MonoisotopicMass(formula = ListFormula("H2O"))
        + MonoisotopicMass(formula = ListFormula("HCOOH")), 
        paste0("loss HCOOH & '", sn$sn, "-H2O' -> MGDG")),
  cbind(MonoisotopicMass(formula = ListFormula("HCOOHCH2")), 
        "loss HCOOH & CH2 -> Lyso PC / PC"),
  c(MonoisotopicMass(formula = ListFormula("HCOOHC7H13NO2")), "loss HCOOH & C7H13NO2 -> Carnitine"),
  cbind(mass2mz(sn$mass, "[M-H]-"), paste0("[", sn$sn, "-H]- -> PA / mPA / PG / PE / Lyso PC / PI / MGDG"))
))
colnames(mzdif.neg) <- c("dif", "add")
mzdif.neg$dif <- as.numeric(mzdif.neg$dif)

## fragments ----
frags <- data.frame(table(round(unlist(mz(ms2)), 3)))
frags <- frags[frags$Freq > 10, ]
frags$Var1 <- as.numeric(as.character(frags$Var1))
mzdif <- get(paste0("mzdif.", tolower(pol)))
mzdif$add <- gsub(" ->.*", "", mzdif$add)
frags$frag <- frags$Var1
for(i in seq(nrow(frags))){
  idx <- which(abs(frags$Var1[i] - mzdif$dif) < 0.01)
  if(length(idx) > 0){
    frags$frag[i] <- paste(mzdif$add[idx], collapse = "; ")
  }
}
dt_f <- data.frame(matrix(ncol = nrow(frags), nrow = length(ms2)))
colnames(dt_f) <- paste0("i_", frags$Var1, " - ", frags$frag)
tmp <- paste(ms2$name, ms2$adduct)
for(i in seq(10)){
  idx <- which(duplicated(tmp))
  if(length(idx) > 0){
    tmp[idx] <- paste0(tmp[idx], letters[i])
  }
}
table(duplicated(tmp))
rownames(dt_f) <- tmp
for(i in seq(nrow(frags))){
  dt_f[,i] <- ifelse(containsMz(ms2, frags$Var1[i], 
                                tolerance = 0.01), 1, 0)
}
idx <- c()
for(i in seq(ncol(dt_f)-1)){
  if(identical(dt_f[,i], dt_f[,i+1]) &
     length(unlist(matchWithPpm(frags$Var1[i], frags$Var1[i+1], ppm = 10))) > 0){
    idx <- c(idx, i)
  }
}
dt_f <- dt_f[,-idx]

## neutral losses ----
ms2_nl <- applyProcessing(ms2)
mz(ms2_nl@backend) <- mz(ms2_nl) - precursorMz(ms2_nl)
nl <- data.frame(table(round(unlist(mz(ms2_nl)), 3)))
nl <- nl[nl$Freq > 5, ]
nl$Var1 <- as.numeric(as.character(nl$Var1))
nl$nl <- nl$Var1
for(i in seq(nrow(nl))){
  idx <- which(abs(abs(nl$Var1[i]) - mzdif$dif) < 0.01)
  if(length(idx) > 0){
    nl$nl[i] <- paste(mzdif$add[idx], collapse = "; ")
  }
}
dt_nl <- data.frame(matrix(ncol = nrow(nl), nrow = length(ms2)))
colnames(dt_nl) <- paste0("l_", -as.numeric(nl$Var1), " - ", nl$nl)
tmp <- paste(ms2$name, ms2$adduct)
for(i in seq(10)){
  idx <- which(duplicated(tmp))
  if(length(idx) > 0){
    tmp[idx] <- paste0(tmp[idx], letters[i])
  }
}
table(duplicated(tmp))
rownames(dt_nl) <- tmp
for(i in seq(nrow(nl))){
  dt_nl[,i] <- ifelse(containsMz(ms2_nl, nl$Var1[i], 
                                 tolerance = 0.01), 1, 0)
}
idx <- c()
for(i in seq(ncol(dt_nl)-1)){
  if(identical(dt_nl[,i], dt_nl[,i+1]) &
     abs(nl$Var1[i] - nl$Var1[i+1]) < 0.01){
    idx <- c(idx, i)
  }
}
dt_nl <- dt_nl[,-idx]

dt <- cbind(dt_f, dt_nl)
dt <- dt[rowSums(dt) > 0,]
dt <- dt[, colSums(dt) > 0]

## network ----
dtx <- dt[,dt["871.5714_19.04",] > 0]
dtx <- dtx[rowSums(dtx) > 2,]
dtx <- dtx[, colSums(dtx) > 2]

net <- graph_from_incidence_matrix(as.matrix(dtx))
V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
V(net)$color[grep("i", names(V(net)))] <- "gold"
V(net)$color[grep("871.5714_19.04", names(V(net)))] <- "black"
V(net)$shape <- c("square", "circle")[V(net)$type+1]
tkplot(net)
plot(net, vertex.frame.color = "White", vertex.label.color = "black")
