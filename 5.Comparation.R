setwd("E:/Rworkspace/ECMO_revise")
library(plyr)
library(timeROC)
library(pROC)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(readxl)
################
# model load   #
################
base_model <- readRDS("./Results/PEEI_model.rds")
final_model<- readRDS("./Results/PEES_model.rds")
summary(base_model)
summary(final_model)

test_imp_data=read.csv("./Results/4.TestSet.csv",check.names = F)
bm_df=read.csv("./Results/2.unicox_data_withoutFill_0.1.csv")
##########################
#  compare to 休克评分   #
##########################
mypal<-pal_lancet("lanonc",alpha=0.6)(9)

colnames(test_imp_data)
tdf=test_imp_data[,c(1:9,13)]
set.seed(699)
death_random=sample(1:373,224)
alive_random=sample(1:648,389)
t1=bm_df$ECMO前休克分层[which(bm_df$event==1)]
t2=bm_df$ECMO前休克分层[which(bm_df$event==0)]
tdf$scai=c(t1[-death_random],t2[-alive_random])
tdf$scai24=test_imp_data$休克分24H


tdf$scai=as.numeric(factor(tdf$scai, levels = c("B", "C", "D", "E"), ordered = TRUE))
tdf$scai24=as.numeric(factor(tdf$scai24, levels = c("B", "C", "D", "E"), ordered = TRUE))
scai=roc(response = tdf$event, predictor = tdf$scai, direction = "<",ci=T)
scai24=roc(response = tdf$event, predictor = tdf$scai24, direction = "<",ci=T)

PEEI.test=test_imp_data$base_risk_score
PEEI.auc=roc(tdf$event,PEEI.test,ci=T)

PEES.test=test_imp_data$risk_score
PEES.auc=roc(tdf$event,PEES.test,ci=T)


plot(PEEI.auc, col = mypal[1], main = "ROC Curves for Multiple Models")
plot(PEES.auc, col = mypal[3], add = TRUE)
plot(scai, col = mypal[5], add = TRUE)
plot(scai24, col = mypal[7], add = TRUE)

# 添加图例
legend("bottomright", legend = c("PEEI Model: AUC=0.730 [0.680,0.780]",
                                 "PEES Model: AUC=0.743 [0.693,0.792]",
                                 "SCAI before ECMO: AUC=0.564 [0.503,0.625]",
                                 "SCAI at 24h ECMO: AUC=0.626 [0.573,0.680]"),
       col = mypal[c(1,3,5,7)], lty = 1)

#DeLong's test专门用于比较两个ROC曲线的AUC差异计算
# DeLong检验 + 多重比较
models <- list(
  PEEI = PEEI.auc,
  PEES = PEES.auc, 
  PECSOS_Pre = scai,
  PECSOS_24h = scai24
)

# 生成所有两两组合
model_names <- names(models)
pairs <- combn(model_names, 2, simplify = FALSE)

# 批量进行DeLong检验
delong_results <- lapply(pairs, function(pair) {
  test <- roc.test(models[[pair[1]]], models[[pair[2]]], method = "delong")
  data.frame(
    Model1 = pair[1],
    Model2 = pair[2],
    AUC1 = round(models[[pair[1]]]$auc, 3),
    AUC2 = round(models[[pair[2]]]$auc, 3),
    P_Value = round(test$p.value, 4)
  )
})

# 合并结果
delong_df <- do.call(rbind, delong_results)

# 多重比较校正
delong_df$P_Bonferroni <- p.adjust(delong_df$P_Value, method = "bonferroni")
delong_df$P_FDR <- p.adjust(delong_df$P_Value, method = "fdr")

cat("=== 完整的DeLong检验结果 (含多重比较校正) ===\n")
print(delong_df)
write.csv(delong_df,"./Results/5.modelCompare_test.csv",row.names = F)

##################
#  validation    #
#################
vali.df=read_xlsx("./Data/2024-2025 清洗模型数据.xlsx")
colnames(vali.df)
vali=data.frame(event=vali.df$辅助最终结局)
vali$event[which(vali$event=="2")]=0 #2=survival 1=dead

library(lubridate)
di=ymd_hm(vali.df$出院时间)-ymd_hm(vali.df$ecmo开始时间)
vali$Time=as.double(di,units="days")
colnames(tdf)
vali$安装前最差PH=vali.df$安装前最差PH
vali$LAC数值=vali.df$安装前最差Lac
#vali.df中 1=否 2=是
vali$ECMO前心脏骤停=0
vali$ECMO前心脏骤停[which(vali.df$ECMO前心脏骤停==2)]=1
vali$心肌炎=0
vali$心肌炎[which(vali.df$入院诊断=="暴发性心肌炎")]=1
vali$出血并发症=0
vali$出血并发症[which(vali.df$出血并发症==2)]=1
vali$肾脏并发症=0
vali$肾脏并发症[which(vali.df$肾并发症==2)]=1
vali$心脏并发症=0
vali$心脏并发症[which(vali.df$心脏并发症==2)]=1

vali$ECMO24h时LAC=vali.df$`24h时-Lac`
vali$scai=vali.df$安装前SCAI分级
vali$scai24=vali.df$`24h-SCAI分级`
vali$scai=as.numeric(factor(vali$scai, levels = c("B", "C", "D", "E"), ordered = TRUE))
vali$scai24=as.numeric(factor(vali$scai24, levels = c("B", "C", "D", "E"), ordered = TRUE))
summary(vali)
write.csv(vali,"./Results/5.validata.csv",row.names = F)

library(autoReg)
library(rrtable)
vali_table=read.csv("./Results/5.validata.csv")

vali_table[,c(5:9)]=as.data.frame(lapply(vali_table[,c(5:9)],as.factor)) 
vali_table$scai=as.factor(vali.df$安装前SCAI分级) 
vali_table$scai24=as.factor(vali.df$`24h-SCAI分级`) 
summary(vali_table)
ft <- gaze(event~.,data=vali_table,digits = 3,method=4) %>%myft()
table2docx(ft,target="./Results/5.Table_vali_iqr.docx")


#### Fill with mean value
vali_mean=read.csv("./Results/5.validata.csv")
vali_mean$安装前最差PH[which(is.na(vali_mean$安装前最差PH))]=7.278
vali_mean$LAC数值[which(is.na(vali_mean$LAC数值))]=6.916
vali_mean$ECMO24h时LAC[which(is.na(vali_mean$ECMO24h时LAC))]=3.443
summary(vali_mean)
scai.vali1=roc(response = vali_mean$event, predictor = vali_mean$scai, direction = "<",ci=T)
scai24.vali1=roc(response = vali_mean$event, predictor = vali_mean$scai24, direction = "<",ci=T)

PEEI.vali1=predict(base_model,newdata = vali_mean,type="risk")
PEEI.auc1=roc(vali_mean$event,PEEI.vali1,ci=T)

PEES.vali1=predict(final_model,newdata = vali_mean,type="risk")
PEES.auc1=roc(vali_mean$event,PEES.vali1,ci=T)


plot(PEEI.auc1, col = mypal[1], main = "ROC Curves for Multiple Models in validation set")
plot(PEES.auc1, col = mypal[3], add = TRUE)
plot(scai.vali1, col = mypal[5], add = TRUE)
plot(scai24.vali1, col = mypal[7], add = TRUE)

# 添加图例
legend("bottomright", legend = c("PEEI Model: AUC=0.775 [0.715,0.834]",
                                 "PEES Model: AUC=0.785 [0.726,0.843]",
                                 "SCAI before ECMO: AUC=0.546 [0.461,0.631]",
                                 "SCAI at 24h ECMO: AUC=0.663 [0.567,0.758]"),
       col = mypal[c(1,3,5,7)], lty = 1)


# DeLong检验 + 多重比较
models <- list(
  PEEI = PEEI.auc1,
  PEES = PEES.auc1, 
  PECSOS_Pre = scai.vali1,
  PECSOS_24h = scai24.vali1
)

# 生成所有两两组合
model_names <- names(models)
pairs <- combn(model_names, 2, simplify = FALSE)

# 批量进行DeLong检验
delong_results <- lapply(pairs, function(pair) {
  test <- roc.test(models[[pair[1]]], models[[pair[2]]], method = "delong")
  data.frame(
    Model1 = pair[1],
    Model2 = pair[2],
    AUC1 = round(models[[pair[1]]]$auc, 3),
    AUC2 = round(models[[pair[2]]]$auc, 3),
    P_Value = round(test$p.value, 4)
  )
})

# 合并结果
delong_df <- do.call(rbind, delong_results)

# 多重比较校正
delong_df$P_Bonferroni <- p.adjust(delong_df$P_Value, method = "bonferroni")
delong_df$P_FDR <- p.adjust(delong_df$P_Value, method = "fdr")

cat("=== 完整的DeLong检验结果 (含多重比较校正) ===\n")
print(delong_df)
write.csv(delong_df,"./Results/5.modelCompare_vali.csv",row.names = F)

