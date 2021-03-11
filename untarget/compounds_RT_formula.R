library(readxl)

z.file <- "C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx"
target <- rbind(read_xlsx(z.file, sheet = "POS"),
                read_xlsx(z.file, sheet = "NEG")) 
colnames(target) <- c("RT", "name")

untarg <- read.csv("compounds.csv")
untarg <- subset(untarg, select = c("name", "formula", "type", "class"))
untarg$name <- gsub("Salt)", "", untarg$name)
untarg$name <- gsub("\\s+$", "", untarg$name)
untarg$name <- gsub(" \\(Na", "", untarg$name)
untarg$name <- gsub(" \\(NH4", "", untarg$name)
untarg$name <- gsub("-H2O", "", untarg$name)
untarg$name <- gsub("enated", "", untarg$name)

data <- merge(target, untarg, by = "name")
write.csv(data, "compounds_RT_formula.csv", row.names = FALSE)


data <- subset(data, class == "IS")
untarg <- subset(untarg, class == "IS")
untarg$name[!untarg$name %in% data$name]


