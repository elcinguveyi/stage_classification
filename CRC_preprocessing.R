library(curatedCRCData)
library(dplyr)

data(GSE12945_eset)
GSE12945 <- as.data.frame(exprs(GSE12945_eset))
GSE12945_pd <- pData(GSE12945_eset)

data(GSE14333_eset)
GSE14333 <- as.data.frame(exprs(GSE14333_eset))
GSE14333_pd <- pData(GSE14333_eset)

data(GSE17536_eset)
GSE17536 <- as.data.frame(exprs(GSE17536_eset))
GSE17536_pd <- pData(GSE17536_eset)

data(GSE17537_eset)
GSE17537 <- as.data.frame(exprs(GSE17537_eset))
GSE17537_pd <- pData(GSE17537_eset)

data(GSE33113_eset)
GSE33113 <- as.data.frame(exprs(GSE33113_eset))
GSE33113_pd <- pData(GSE33113_eset)

data(GSE39582_eset)
GSE39582 <- as.data.frame(exprs(GSE39582_eset))
GSE39582_pd <- pData(GSE39582_eset)

data(GSE24549.GPL5175_eset)
GSE24549.GPL5175 <- as.data.frame(exprs(GSE24549.GPL5175_eset))
GSE24549.GPL5175_pd <- pData(GSE24549.GPL5175_eset)

data(GSE24550.GPL5175_eset)
GSE24550.GPL5175 <- as.data.frame(exprs(GSE24550.GPL5175_eset))
GSE24550.GPL5175_pd <- pData(GSE24550.GPL5175_eset)


common_genes <- intersect(intersect(intersect(intersect(intersect(intersect(intersect(rownames(GSE12945),rownames(GSE14333)),rownames(GSE17536)),rownames(GSE17537)),rownames(GSE33113)),rownames(GSE39582)),rownames(GSE24549.GPL5175)),rownames(GSE24550.GPL5175))

GSE12945 <- GSE12945[rownames(GSE12945) %in% common_genes, ] 
GSE14333 <- GSE14333[rownames(GSE14333) %in% common_genes, ] 
GSE17536 <- GSE17536[rownames(GSE17536) %in% common_genes, ]
GSE17537 <- GSE17537[rownames(GSE17537) %in% common_genes, ]
GSE33113 <- GSE33113[rownames(GSE33113) %in% common_genes, ] 
GSE39582 <- GSE39582[rownames(GSE39582) %in% common_genes, ]
GSE24549.GPL5175 <- GSE24549.GPL5175[rownames(GSE24549.GPL5175) %in% common_genes, ]
GSE24550.GPL5175 <- GSE24550.GPL5175[rownames(GSE24550.GPL5175) %in% common_genes, ] 

integrated <-as.data.frame(rbind(t(GSE12945),t(GSE14333)))
integrated <-as.data.frame(rbind(integrated,t(GSE17536)))
integrated <-as.data.frame(rbind(integrated,t(GSE17537)))
integrated <-as.data.frame(rbind(integrated,t(GSE33113)))
integrated <-as.data.frame(rbind(integrated,t(GSE39582)))
integrated <-as.data.frame(rbind(integrated,t(GSE24549.GPL5175)))
integrated <-as.data.frame(rbind(integrated,t(GSE24550.GPL5175)))

# Add Dataset Name Column
GSE12945_pd["dt_name"] <- "GSE12945"
GSE14333_pd["dt_name"] <- "GSE14333"
GSE17536_pd["dt_name"] <- "GSE17536"
GSE17537_pd["dt_name"] <- "GSE17537"
GSE33113_pd["dt_name"] <- "GSE33113"
GSE39582_pd["dt_name"] <- "GSE39582"
GSE24549.GPL5175_pd["dt_name"] <- "GSE24549.GPL5175"
GSE24550.GPL5175_pd["dt_name"] <- "GSE24550.GPL5175"

integrated_pd <- GSE12945_pd[,c(3,29,31,32,60)]
integrated_pd <- rbind(integrated_pd,GSE14333_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE17536_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE17537_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE33113_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE39582_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE24549.GPL5175_pd[,c(3,29,31,32,60)])
integrated_pd <- rbind(integrated_pd,GSE24550.GPL5175_pd[,c(3,29,31,32,60)])


# Remove Normal Samples, NA and stage-0 samples
integrated <- integrated[which(integrated_pd$sample_type == "tumor"),]
integrated_pd <- integrated_pd[which(integrated_pd$sample_type == "tumor"),]
integrated <- integrated[which(is.na(integrated_pd$stageall)==FALSE),]
integrated_pd <- integrated_pd[which(is.na(integrated_pd$stageall)==FALSE),]
integrated <- integrated[-c(which(integrated_pd$stageall==0)),]
integrated_pd <- integrated_pd[-c(which(integrated_pd$stageall==0)),]

write.csv(integrated,file = "8datasets-Data/integrated.csv")
write.csv(integrated_pd,file = "8datasets-Data/integrated_pd.csv")

# BATCH EFFECT REMOVAL

library(ggfortify)
library(limma)

# PCA - Before
pca_before <- prcomp(integrated, scale. = TRUE)

integrated["dt_name"] <- integrated_pd[,5]
svg("8datasets-Res/pca-before.svg")

autoplot(pca_before,data=integrated,  colour = 'dt_name')
dev.off()

integrated <- integrated[,-c(6508)] #remove dataset name column
integrated.mat<-as.matrix.data.frame(t(integrated))
mode(integrated.mat)='numeric'

integrated.corrected <- removeBatchEffect(integrated.mat, integrated_pd[,5])
integrated.corrected<-as.data.frame(t(integrated.corrected))
colnames(integrated.corrected)<-colnames(integrated)
rownames(integrated.corrected)<-rownames(integrated)

#PCA - After
pca_after <- prcomp(integrated.corrected, scale. = TRUE)

integrated.corrected["dt_name"] <- integrated_pd[,5]
svg("8datasets-Res/pca-after.svg")

autoplot(pca_after, data=integrated.corrected, colour = 'dt_name')
dev.off()

integrated.corrected <- integrated.corrected[,-c(6508)] #remove dataset name column

par(mar = rep(2, 4))
png("8datasets-Res/boxplot-before.png",width = 10000,height = 5000, res=600)
boxplot(t(integrated),main="Original")
dev.off()
png("8datasets-Res/boxplot-after.png",width = 10000,height = 5000, res=600)
boxplot(t(integrated.corrected),main="Batch Corrected")
dev.off()

write.csv(integrated.corrected,file = "8datasets-Data/integrated_corrected.csv")
###


