untarg <- read.csv("compounds.csv")
target <- read.csv("../target/data/target_detected_RT_MRM.csv")
target$ID <- gsub("-CH3", "", target$ID)
target$ID <- gsub("_Na", "", target$ID)
for(i in 1:nrow(untarg)){
  if(untarg$RT_checked[i] == 0){
    if(untarg$name[i] %in% target$ID){
      untarg$RT[i] <- round(target$RT[which(target$ID == untarg$name[i])]*60) + 12
    }
  }
}
write.csv(untarg, "compounds.csv", row.names = F)


untarg$RT_target <- NA
for(i in seq(nrow(untarg))){
  if(untarg$name[i] %in% target$ID){
    untarg$RT_target[i] <- round(target$RT[which(target$ID == untarg$name[i])]*60)
  }
}
untarg$RT_delta <- untarg$RT - untarg$RT_target
mean(untarg$RT_delta, na.rm = T)
