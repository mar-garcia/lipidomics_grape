library(readxl)

z.file <- "C:/Users/garciaalom/Google Drive/projectes/lipidomics_shared/new_rt.xlsx"
target <- rbind(read_xlsx(z.file, sheet = "POS"),
                read_xlsx(z.file, sheet = "NEG")) 
colnames(target) <- c("RT", "ID")

form <- rbind(read_xlsx(z.file, sheet = "FORMULE"),
              read_xlsx(z.file, sheet = "IS"))

data <- merge(target, form, by = "ID")
write.csv(data, "compounds_RT_formula.csv", row.names = FALSE)


data <- subset(data, class == "IS")
form <- subset(form, class == "IS")
form$ID[!form$ID %in% data$ID]


