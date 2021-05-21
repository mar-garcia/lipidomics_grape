library(readxl)

z.file <- "C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx"
target <- rbind(read_xlsx(z.file, sheet = "POS"),
                read_xlsx(z.file, sheet = "NEG")) 
colnames(target) <- c("RT", "ID")
target$ID <- gsub("_Na", "", target$ID)
target$ID <- gsub("-H2O", "", target$ID)
target$ID <- gsub("-COO", "", target$ID)
target$ID[target$ID == "lact_dihydroCER_24:0"] <- "lactCER_24:0"

form <- rbind(read_xlsx(z.file, sheet = "FORMULE"),
              read_xlsx(z.file, sheet = "IS"))
form$ID <- gsub("_Na", "", form$ID)

data <- merge(target, form, by = "ID")
write.csv(data, "compounds_RT_formula.csv", row.names = FALSE)


data <- subset(data, class == "IS")
form <- subset(form, class == "IS")
form$ID[!form$ID %in% data$ID]


