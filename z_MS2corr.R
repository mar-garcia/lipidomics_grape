library(MsFeatures)
library(MetaboCoreUtils)
library(tidyverse)
library(OrgMassSpecR)
myf <- function(x){
  SpectrumSimilarity(
    x[[1]], x[[2]], 
    print.graphic = FALSE, t = 0.05, b = 20)
}
load("data/RData.RData")
load("output/pks.RData")

for(s in c("tissues", "maturation")){
  for(p in c("POS", "NEG")){
    pks <- get(paste("pks", s, p, sep = "_"))
    pks$FT <- rownames(pks)
    assign(paste("pks", s, p, sep = "_"), pks)
  }
}
pks <- rbind(pks_maturation_NEG, pks_maturation_POS, 
             pks_tissues_NEG, pks_tissues_POS)
pks <- pks[!is.na(pks$score_rt),]
pks <- pks[!grepl("-", pks$target_compound_id),]
pks <- pks[!grepl("2M", pks$target_ion_adduct),]
pks <- pks[pks$n_MS2 > 1, ]
pks <- pks[order(pks$n_MS2), ]
pks$sim <- NA
for(i in seq(nrow(pks))){
  dda_xdata <- get(paste("dda_xdata", pks$dataset[i], sep = "_"))
  dda_spectra <- get(paste("dda_spectra", pks$dataset[i], sep = "_"))
  ft <- get(paste("ft", pks$dataset[i], sep = "_"))
  i_id <- rownames(
    chromPeaks(dda_xdata)[ft[pks$FT[i], "peakidx"][[1]],])
  i_spectra <- dda_spectra[dda_spectra$peak_id %in% i_id]
  mypairs <- t(combn(seq(length(i_spectra)), 2))
  pippo <- mypairs %>%
    as_tibble() %>% 
    mutate(mylist = map2(V1,V2, function(x,y) c(x,y))) %>% 
    pull(mylist)
  pluto <- lapply(pippo, function(x) {
    list(cbind( "mz" = unlist(mz(i_spectra[x[1]])), 
                "intensity" = unlist(intensity(i_spectra[x[1]]))),
         cbind( "mz" = unlist(mz(i_spectra[x[2]])), 
                "intensity" = unlist(intensity(i_spectra[x[2]]))))
  }) 
  pks$sim[i] <- paste(sprintf(
    "%.3f", range(unlist(lapply(pluto, myf)))), collapse = "-")
  write.csv(pks, "output/MS2sim.csv")
}
