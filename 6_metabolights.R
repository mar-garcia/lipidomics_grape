load("output/identifications.RData")
is <- 1:grep("d7", cmps$name)[length(grep("d7", cmps$name))]
cmps <- cmps[-is,]
cmpsdb <- read.csv("data/compounds.csv")

for(p in c("NEG", "POS")){
  if(p == "POS"){
    tmp_cmps <- cmps[cmps$name %in% pks$target_name[grep(p, pks$dataset)] |
                       grepl("NH4", cmps$tissues) |
                       grepl("NH4", cmps$maturation),]
    tmp_cmps$charge <- 1
  } else if(p == "NEG"){
    tmp_cmps <- cmps[cmps$name %in% pks$target_name[grep(p, pks$dataset)],]
    tmp_cmps$charge <- -1
  }
  tmp_cmps <- tmp_cmps[c("name", "formula", "RT", "charge", "class")]
  tmp_cmps$mz <- NA
  tmp_cmps$adduct <- NA
  tmp_pks <- pks[grep(p, pks$dataset),]
  
  for(i in seq(nrow(tmp_cmps))){
    i.add <- unique(tmp_pks$target_ion_adduct[
      tmp_pks$target_name == tmp_cmps$name[i]])
    if((length(i.add) == 0) & (tmp_cmps$class[i] == "TG")){
      i.add <- unique(tmp_pks$target_ion_adduct[
        tmp_pks$target_name == cmpsdb$compound[
          cmpsdb$comments == tmp_cmps$name[i]]])
    }
    if(length(i.add) > 1){
      i.add <- i.add[!grepl("2M", i.add)]
    }
    if((length(i.add) > 1) & (p == "NEG") & (tmp_cmps$class[i] == "Cer;O4")){
      i.add <- i.add[grepl("CHO2", i.add)]
    }else if((length(i.add) > 1) & (p == "NEG") & 
             (tmp_cmps$class[i] == "DGMG")){
      i.add <- i.add[!grepl("CHO2", i.add)]
    }else if((length(i.add) > 1) & (p == "NEG") & 
             (grepl("Malic|catechin", tmp_cmps$name[i]))){
      i.add <- "[M-H]-"
    }
    if(length(i.add) != 1){
      print(i)
    }
    tmp_pks2 <- tmp_pks[(tmp_pks$target_name == tmp_cmps$name[i]) & 
                          (tmp_pks$target_ion_adduct == i.add), ]
    if(nrow(tmp_pks2) == 0){
      tmp_pks2 <- tmp_pks[
        (tmp_pks$target_name == cmpsdb$compound[
          cmpsdb$comments == tmp_cmps$name[i]]) & 
          (tmp_pks$target_ion_adduct == i.add), ]
    }
    idx <- which.min(abs(tmp_pks2$ppm_error))
    tmp_cmps$mz[i] <- tmp_pks2$mzmed[idx]
    tmp_cmps$adduct[i] <- tmp_pks2$target_ion_adduct[idx]
  }
  tmp_cmps$mz <- round(tmp_cmps$mz, 4)
  write.csv(tmp_cmps, paste0("output/MTBLS_", p, ".csv"), row.names = FALSE)
}

