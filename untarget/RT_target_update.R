setwd("~/GitHub/lipidomics/untarget")
target <- read.csv("../target/data/target_detected_RT_MRM.csv")
untarg <- read.csv("data/compounds.csv")
for(i in 1:nrow(untarg)){
  if(untarg$RT_checked[i] == 0){
    if(untarg$name[i] %in% target$ID){
      untarg$RT[i] <- round(target$RT[which(target$ID == untarg$name[i])]*60)
    }
  }
}
write.csv(untarg, "data/compounds.csv")
