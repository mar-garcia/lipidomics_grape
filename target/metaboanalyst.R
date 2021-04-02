library(MetaboAnalystR)
mSet<-InitDataObjects("conc", "ts", FALSE)
mSet<-SetDesignType(mSet, "time0")
mSet<-Read.TextData(mSet, "test.csv", "rowts", "disc")
#mSet<-SanityCheckData(mSet)
#mSet<-ContainMissing(mSet)
#mSet<-ReplaceMin(mSet)
#mSet<-PreparePrenormData(mSet)
#mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
#mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
#mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-ANOVA2.Anal(mSet, 0.05, "fdr", "time0", 1, 1)
tmp <- data.frame(mSet$analSet$aov2$sig.mat)
tmp <- tmp[order(tmp$Adjusted.P.val),]
head(tmp)
PlotCmpdSummary(mSetObj = mSet, cmpdNm = "4.3148/653", #0, 
                      format="png", dpi=72, width=NA)
mSet <- performMB(mSet, 10)
tmp <- data.frame(mSet$analSet$MB$stats)
tmp$id <- rownames(tmp)
tmp <- tmp[order(tmp$Hotelling.T2, decreasing = T),]
head(tmp)
PlotMBTimeProfile(mSet, "1.8927/853", "png", 72, width = NA)



########################################

library(MetaboAnalystR)
mSet<-InitDataObjects("conc", "ts", FALSE)
mSet<-SetDesignType(mSet, "time0")
mSet<-Read.TextData(mSet, "data/dt_metaboanalyst.csv", "rowts", "disc");
#mSet<-SanityCheckData(mSet)
#mSet<-ContainMissing(mSet)
#mSet<-ReplaceMin(mSet)
#mSet<-PreparePrenormData(mSet)
#mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
#mSet<-PlotNormSummary(mSet, "norm_1_", "png", 72, width=NA)
#mSet<-PlotSampleNormSummary(mSet, "snorm_1_", "png", 72, width=NA)
mSet<-ANOVA2.Anal(mSet, 2, "fdr", "time0", 1, 1)
mSet <- performMB(mSet, ncol(mSet$dataSet$orig))

tmp <- data.frame(mSet$analSet$aov2$sig.mat)
tmp2 <- data.frame(mSet$analSet$MB$stats)


write.csv(tmp, "pval_metaboanalyst.csv")
