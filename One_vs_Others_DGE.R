library(limma)
library(tidyr)
library(Biobase)


df_main <- read.csv("stage_vs_others_dge.csv")
df_surv <- read.csv("integrated_pd.csv")
df <- df_main
df$stage=df_surv$stageall


df$stage[df$stage==0] = NA  #stage "0" removing.
df = data.frame(df)
df = df[!is.na(df$stage),]


rownames(df) <- df$X  #first column was named as X. so we convert this to rownames and then delete this column.
df$X = NULL

df$stage[df$stage != 1] = 2
p = df$stage   #phenodata 
df$stage = NULL
x = df

class(x)
x = as.matrix(x)  #assaydata should be in matrix form.
dim(x)
x = t(x)

class(p)
p = data.frame(p)
dim(p)
colnames(p) 
rownames(p) = rownames(df)
rownames(p)
head(p)


colnames(p)[which(names(p) == "p")] <- "Label"
p$Label[p$Label==1] = 'Stage1' 
p$Label[p$Label==2] = 'StageOther'
#Bu kod, stage1 vs diðerleri olacak þekilde çalýþmaktadýr. Diðer stageler için kodu çalýþtýrmak istediðimizde,
#örneðin stage 2 için yapacaksak, p$Label[p$Label==2] = 'Stage2',  p$Label[p$Label==3] = 'StageOther' þeklinde düzenleriz.

eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p))
table(pData(eset)[, "Label"])

design <- model.matrix(~0 + Label, data = pData(eset))                      
head(design, 4)
colSums(design)

cm <- makeContrasts(S2vOther = LabelStage1 - LabelStageOther,
                    levels = design)


# Fit coefficients
fit <- lmFit(eset, design)
# Fit contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)
# Calculate t-statistics
fit2 <- eBayes(fit2)
# Summarize results
results <- decideTests(fit2) #default: p.value=0.05,lfc=0 
summary(results)


r = data.frame(results)
length(rownames(r))

genes = vector()
for (i in 1:nrow(r)){
  if(r[i,1] != 0  )
    genes = c(genes, rownames(r)[i] )
}
gene_list_1 = unique(genes)

#bu þekilde her bir stage için ayrý ayrý gen listleri elde edildikten sonra listeler birleþtirilir:
all_genes=vector()
all_genes=c(gene_list_1, gene_list_2, gene_list_3, gene_list_4)
all_genes=unique(all_genes)  
write(allgenes, "allGenes.txt")


#var olan verisetinde sadece bu genler kalacak þekilde düzenleme yapýp kaydetmek istersek:
gene_list = all_genes
gene_list = as.data.frame(gene_list) 
new_dataset = subset(df, select=gene_list$gene_list)
write.csv(new_dataset, file = "new_dataset.csv")


