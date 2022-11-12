library("survival")
library("survminer")
library("readxl")
library(survival)
library(survminer)
library(matrixStats)
library(readxl)
library(curatedCRCData)
library(dplyr)

HR_anly <- function(hr_data,label,type,stage){
  #Univarita Analysis with multiple covariates *******
  path_hr = paste0("8datasets-Res/HR/")
  
  hr_data <- hr_data[which(label == stage),]
  
  covariates <-  colnames(hr_data[,3:ncol(hr_data)])#HR analizini yapmak istedigin genleri covariate olarak belirle
  
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time = as.numeric(hr_data$days_to_recurrence_or_death), event = as.numeric(hr_data$dfs_status))~', x)))
  
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = hr_data)})
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           res<-c(beta, HR,HR.confint.lower,HR.confint.upper, wald.test, p.value)
                           names(res)<-c("beta", "HR","Low","High", "wald.test", 
                                         "p.value")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res)# HR Analiz sonuclari
  res
  
  write.csv(res,file=paste0(path_hr,type,"_hr", ".csv"))
  # library(forestplot)
  # png(paste0(path_hr,"_",type,"_hr", ".png"), width = 20, height = 20, pointsize=12)
  # forestplot( rownames(res), as.numeric(res$HR),as.numeric(res$Low), as.numeric(res$High), xlog=T, boxsize = 0.2,ci.vertices = T) #grafik
  # dev.off()
}


s3_downGenes <- c("ARSE","CHN2", "EPHB2", "ACOX2")
s4_downGenes <- c("CD177")
s3_upGenes <- c("SERPINB5","SCEL", "NTM", "CDKN2B","AKAP12","AHNAK2")
s4_upGenes <- c("PLN", "SPP1", "MGP", "PEG10","AKAP12")


load("8datasets-Data/integrated_corrected.rda")
df <- integrated_corrected


load("8datasets-Data/integrated_pd.rda") #survival verisini oku
surv_data <- integrated_pd[,c(2:4)] #DFS sutunlarini aldim
surv_data[which(surv_data$dfs_status == "deceased_or_recurrence"),2] = 1
surv_data[which(surv_data$dfs_status == "living_norecurrence"),2] = 0

surv_data <- surv_data[which(is.na(surv_data[,2])==FALSE),]
surv_data <- surv_data[which(is.na(surv_data[,3])==FALSE),]

reorder_idx <- match(rownames(surv_data),rownames(df)) # Saving indices for how to reorder `second` to match `first`
df <- df[reorder_idx,] 

s3_downExp <- df[, which((colnames(df) %in% s3_downGenes)==TRUE)]
s4_downExp <- df[, which((colnames(df) %in% s4_downGenes)==TRUE)]
s3_upExp <- df[, which((colnames(df) %in% s3_upGenes)==TRUE)]
s4_upExp <- df[, which((colnames(df) %in% s4_upGenes)==TRUE)]

#tek gen için covariates satiri düzeltilmeli

HR_anly(cbind(surv_data[,c(2,3)],s4_downExp),surv_data[,1],"Stage4_down",4)










