rtimes <- c()
ppms <- c()
annots <- c()

for(i in 1:nrow(features)){
  rtr <- c(features$rtmed[i], features$rtmed[i]) + 10 * c(-1, 1)
  cmp.i <- ions_long[unlist(
    matchWithPpm(features$mzmed[i], ions_long$mz, ppm = 10)), ]
  cmp.i <- cmp.i[cmp.i$RT > rtr[1] & cmp.i$RT < rtr[2], ]
  rtimes <- c(rtimes, features$rtmed[i] - cmp.i$RT)
  ppms <- c(ppms, ((features$mzmed[i] - cmp.i$mz)/cmp.i$mz)*1e6)
  annots <- c(annots, cmp.i$annotation)
}


plot(rtimes[order(rtimes)])
abline(h = c(-1.4, 2), col = "red") # maturation POS
plot(rtimes[order(rtimes)], ylim = c(-1.3, 1.8)) # maturation NEG

plot(ppms[order(ppms)])
abline(h = c(-8, -0.7), col = "red") # maturation POS
plot(ppms[order(ppms)], ylim = c(-2.3,-0.2)) # maturation NEG


data.frame(table(annots))[order(data.frame(table(annots))[,2]),]
