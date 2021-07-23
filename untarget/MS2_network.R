library(igraph)
library(pheatmap)
load("~/GitHub/lipidomics/untarget/output/MS2_spectra.RData")
p <- "NEG"
sps_ms2 <- get(paste0("sps_ms2_", p))
sps_ms2 <- replaceIntensitiesBelow(sps_ms2, threshold = 10, value = 0)
sps_ms2 <- filterIntensity(sps_ms2, intensity = c(0.1, Inf))

## fragments ----
frags <- data.frame(table(round(unlist(mz(sps_ms2)), 2)))
frags <- frags[frags$Freq >= 5, ]
#frags$Var1 <- as.numeric(as.character(frags$Var1))
#mzdif <- get(paste0("mzdif.", tolower(pol)))
#mzdif$add <- gsub(" ->.*", "", mzdif$add)
#frags$frag <- frags$Var1
#for(i in seq(nrow(frags))){
#  idx <- which(abs(frags$Var1[i] - mzdif$dif) < 0.01)
#  if(length(idx) > 0){
#    frags$frag[i] <- paste(mzdif$add[idx], collapse = "; ")
#  }
#}
dt_f <- data.frame(matrix(ncol = nrow(frags), nrow = length(sps_ms2)))
colnames(dt_f) <- paste0("i_", frags$Var1)
rownames(dt_f) <- sps_ms2$name
#for(i in seq(10)){
#  idx <- which(duplicated(tmp))
#  if(length(idx) > 0){
#    tmp[idx] <- paste0(tmp[idx], letters[i])
#  }
#}
#table(duplicated(tmp))
#rownames(dt_f) <- tmp
for(i in seq(nrow(frags))){
  dt_f[,i] <- ifelse(containsMz(sps_ms2, as.numeric(as.character(frags$Var1[i])), 
                                tolerance = 0.1), 1, 0)
}
#idx <- c()
#for(i in seq(ncol(dt_f)-1)){
#  if(identical(dt_f[,i], dt_f[,i+1]) &
#     length(unlist(matchWithPpm(frags$Var1[i], frags$Var1[i+1], ppm = 10))) > 0){
#    idx <- c(idx, i)
#  }
#}
#dt_f <- dt_f[,-idx]

## neutral losses ----
ms2_nl <- applyProcessing(sps_ms2)
mz(ms2_nl@backend) <- mz(ms2_nl) - precursorMz(ms2_nl)
nl <- data.frame(table(round(unlist(mz(ms2_nl)), 2)))
nl <- nl[nl$Freq >= 5, ]
#nl$Var1 <- as.numeric(as.character(nl$Var1))
#nl$nl <- nl$Var1
#for(i in seq(nrow(nl))){
#  idx <- which(abs(abs(nl$Var1[i]) - mzdif$dif) < 0.01)
#  if(length(idx) > 0){
#    nl$nl[i] <- paste(mzdif$add[idx], collapse = "; ")
#  }
#}
dt_nl <- data.frame(matrix(ncol = nrow(nl), nrow = length(ms2_nl)))
colnames(dt_nl) <- paste0("l_", -as.numeric(nl$Var1))
rownames(dt_nl) <- ms2_nl$name
#tmp <- paste(ms2$name, ms2$adduct)
#for(i in seq(10)){
#  idx <- which(duplicated(tmp))
#  if(length(idx) > 0){
#    tmp[idx] <- paste0(tmp[idx], letters[i])
#  }
#}
#table(duplicated(tmp))
#rownames(dt_nl) <- tmp
for(i in seq(nrow(nl))){
  dt_nl[,i] <- ifelse(containsMz(ms2_nl, as.numeric(as.character(nl$Var1[i])), 
                                 tolerance = 0.1), 1, 0)
}
#idx <- c()
#for(i in seq(ncol(dt_nl)-1)){
#  if(identical(dt_nl[,i], dt_nl[,i+1]) &
#     abs(nl$Var1[i] - nl$Var1[i+1]) < 0.01){
#    idx <- c(idx, i)
#  }
#}
#dt_nl <- dt_nl[,-idx]

dt <- cbind(dt_f, dt_nl)
dt <- dt[rowSums(dt) > 0,]
dt <- dt[, colSums(dt) > 0]

## network ----
#dtx <- dt[,dt["869.5587_19.01",] > 0]
#dtx <- dtx[rowSums(dtx) > 2,]
#dtx <- dtx[, colSums(dtx) > 2]
dtx <- dt

net <- graph_from_incidence_matrix(as.matrix(dtx))
V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
#V(net)$color[grep("i", names(V(net)))] <- "gold"
#V(net)$color[grep("871.5714_19.04", names(V(net)))] <- "black"
V(net)$shape <- c("square", "circle")[V(net)$type+1]
tkplot(net)
plot(net, vertex.frame.color = "White", vertex.label.color = "black")




res <- compareSpectra(sps_ms2, tolerance = 0.005)
colnames(res) <- rownames(res) <- sps_ms2$name
pheatmap(res, main = "Dotproduct")

res <- compareSpectra(ms2_nl, tolerance = 0.005)
colnames(res) <- rownames(res) <- ms2_nl$name
pheatmap(res, main = "Dotproduct")
