library(survival)
library(survminer)
library(matrixStats)
library(readxl)
library(curatedCRCData)
library(dplyr)

survival_plot <- function(data, surv_data,dataset,label,stage_reguleGenes,type,stage){
  path_surv = paste0("8datasets-Res/survival/",dataset,"/stageVsOther/")
  
  time_ ="days_to_recurrence_or_death" #"days_to_death"#"DFS.Time"#"days_to_recurrence_or_death" 
  ind_ = "dfs_status"#"vital_status"#"DFS"#"dfs_status"
  
  
  splots <- list()
  for( i in 1:20){
    
    filt_genes <- match(stage_reguleGenes[i,1],colnames(data),nomatch = 0) 
    data.filt <- as.data.frame(data[,filt_genes]) 
    data.filt <- as.data.frame(data.filt[which(label == stage),])
    surv_stage <- surv_data[which(label == stage),]
    #surv_stage <- surv_data
    
    xvar <- rowMedians(as.matrix(data.filt))
    time <- as.numeric(surv_stage[,time_])
    censor <- as.numeric(surv_stage[,ind_])
    dat <- cbind.data.frame(time, censor, xvar)
    cat <- dat
    th <- median(cat[,3])
    cat[,3][dat[,3] >= th ] <- "H"
    cat[,3][dat[,3] < th] <- "L"
    time <- cat[,paste("time", collapse = NULL, sep = "")]
    censor <- cat[,paste("censor", collapse = NULL, sep = "")]
    table(cat[,3])#control high-low size
    #bar_plot(cat,label)
    fit <- survfit(Surv(as.numeric(time), as.numeric(factor(censor))) ~ xvar,  data = cat)
    
    splots[[i]] <- ggsurvplot(fit, data = cat,
                              #conf.int = TRUE,
                              surv.median.line = "hv",
                              tables.y.text = FALSE,
                              tables.theme = theme_cleantable(),
                              xlab = "Time in days",
                              pval = T,
                              # pval.method = TRUE,
                              pval.coord = c(0.05, 0.05),
                              risk.table = TRUE,
                              legend = c(0.9, 0.2), #c("bottom"), #
                              legend.labs = c("High","Low"),
                              legend.title = c("Group"),
                              title=paste0(stage_reguleGenes[i,1]),
                              #title=paste0("Stage1&2"),
                              fontsize = 4,
                              ggtheme = theme_bw(base_family = "sans", base_size = 10))
    
  }
  png(paste0(path_surv,dataset,"_",type,"_survival", ".png"), width = 4000, height = 600, pointsize=12)
  p.surv <- arrange_ggsurvplots(splots, print = FALSE,labels = LETTERS[1:10], ncol = 10, nrow = 2)
  
  print(p.surv)
  dev.off()
  
  svg(paste0(path_surv,dataset,"_",type,"_survival", ".svg"), width = 25, height = 5, pointsize=12)
  print(p.surv)
  dev.off()
  
  
  
  
}
dataset = "Integrated"
df <- integrated_corrected

df_label <- as.data.frame(integrated_pd[,2])
rownames(df_label) <- rownames(integrated_pd)
colnames(df_label) <- "Label"

df_survival <- integrated_pd[,c(3,4)]
df_survival[which(df_survival$dfs_status == "deceased_or_recurrence"),1] = 1
df_survival[which(df_survival$dfs_status == "living_norecurrence"),1] = 0

df_survival <- df_survival[which(is.na(df_survival[,1])==FALSE),]
df_survival <- df_survival[which(is.na(df_survival[,2])==FALSE),]

reorder_idx <- match(rownames(df_survival),rownames(df)) # Saving indices for how to reorder `second` to match `first`
df <- df[reorder_idx,] 
df_label <- df_label[reorder_idx,]

type ="down" #down
stage1_reguleGenes <- read_excel(paste0("8datasets-Data/stageVsOther/",type,"_regule_genes.xlsx"), 
                                 sheet = "Stage1")
#stage2_reguleGenes <- read_excel(paste0("8datasets-Data/stageVsOther/",type,"_regule_genes.xlsx"), 
 #                                sheet = "Stage2") 
stage3_reguleGenes <- read_excel(paste0("8datasets-Data/stageVsOther/",type,"_regule_genes.xlsx"), 
                                 sheet = "Stage3")
stage4_reguleGenes <- read_excel(paste0("8datasets-Data/stageVsOther/",type,"_regule_genes.xlsx"), 
                                 sheet = "Stage4")


survival_plot(df,df_survival,dataset,df_label,stage4_reguleGenes,"Stage4_up",4)


