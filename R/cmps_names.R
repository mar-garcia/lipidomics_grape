
cmps <- cdb_RT

cmps_names <- function(cmps){
  cmps <- cmps[!cmps$compound %in% c(
    "FFA(22:0)", "PA(18:1/18:3)", "PC(18:1/18:3)", 
    "PE(18:1/18:1)", "PE(18:1/18:3)"),]
  
  cmps$class <- gsub("DAG", "DG", cmps$class)
  cmps$class <- gsub("TAG", "TG", cmps$class)
  cmps$class <- gsub("CER_hex_dihydroxy", "HexCer;O4", cmps$class)
  cmps$class <- gsub("CER_hex_hydroxy", "HexCer;O3", cmps$class)
  cmps$class <- gsub("CER_hex", "HexCer", cmps$class)
  cmps$class <- gsub("CER_dihydroxy", "Cer;O4", cmps$class)
  cmps$class <- gsub("CER_hydroxy", "Cer;O3", cmps$class)
  cmps$class <- gsub("CER", "Cer", cmps$class)
  cmps$class <- gsub("non_FA", "others", cmps$class)
  cmps$class <- gsub("sitosterol_hex", "sitosterol_ester_hex", cmps$class)
  cmps$class <- gsub("Triterpenoids", "triterpen", cmps$class)
  cmps$class <- gsub("unknown", "unk", cmps$class)
  
  cmps$compound <- gsub("\\(", "\\ ", cmps$compound)
  cmps$compound <- gsub(")", "", cmps$compound)
  cmps$compound <- gsub("d4", "\\(d4)", cmps$compound)
  cmps$compound <- gsub("d7", "\\(d7)", cmps$compound)
  cmps$compound <- gsub("CER", "Cer", cmps$compound)
  cmps$compound <- gsub("Cer_hex", "HexCer", cmps$compound)
  cmps$compound <- gsub("Cer_hydroxy", "Cer;O3", cmps$compound)
  cmps$compound <- gsub("Cer_dihydroxy", "Cer;O4", cmps$compound)
  cmps$compound <- gsub("DAG", "DG", cmps$compound)
  cmps$compound <- gsub("TAG", "TG", cmps$compound)
  cmps$compound[grep("catechin", cmps$compound)] <- "Epicatechin"
  #cmps$compound <- gsub("\\ IV", " \\(IV)", cmps$compound)
  #cmps$compound <- gsub("\\ III", " \\(III)", cmps$compound)
  #cmps$compound <- gsub("\\ II", " \\(II)", cmps$compound)
  #cmps$compound <- gsub("\\ I", " \\(I)", cmps$compound)
  return(cmps)
}
