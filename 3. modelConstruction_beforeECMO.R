setwd("E:/Rworkspace/ECMO_revise")
library(dplyr)
library(survival)
library(survminer)
library(plyr)
library(MASS)
library(pROC)
library(timeROC)
library(rms)
library(ggplot2)
library(ggpubr)
library(readxl)
library(impute)
#########################
#   before ECMO         #
#########################
#Prognosis Evaluation at ECMO Initiation (PEEI) model
#mix imp and mucox
data=read.csv("./Results/2.final_Data.csv")
colnames(data)
summary(data)
data$Time=as.numeric(data$Time)
colnames(data)
# data$LAC数值[which(is.na(data$LAC数值))]=6.933
# data$安装前最差PH[which(is.na(data$安装前最差PH))]=7.282
#分训练测试集6：4
death=data[which(data$event==1),]
alive=data[which(data$event==0),]
set.seed(699)
death_random=sample(1:373,224)
alive_random=sample(1:648,389)

#training set
train=rbind(death[death_random,],alive[alive_random,])
#fill NA with knn in the train set
data_imputed <- as.data.frame(impute.knn(as.matrix(train[, 1:9]))$data)

comparison_results <- data.frame()
for (var in colnames(train)[3:9]) {
  # 原始数据（去除缺失值）
  original <- na.omit(train[[var]])
  imputed <- data_imputed[[var]]
  
  # 计算统计量
  stats <- data.frame(
    Variable = var,
    Original_Mean = round(mean(original), 3),
    Imputed_Mean = round(mean(imputed), 3),
    Original_SD = round(sd(original), 3),
    Imputed_SD = round(sd(imputed), 3),
    Mean_Diff = round(abs(mean(original) - mean(imputed)), 3),
    N_Observed = length(original),
    N_Imputed = sum(is.na(train[[var]]))
  )
  
  comparison_results <- rbind(comparison_results, stats)
}
comparison_results$relative_diff <- comparison_results$Mean_Diff / comparison_results$Original_SD
print(comparison_results)

train<-data_imputed
table(test$event)
#test set
test=rbind(death[-death_random,],alive[-alive_random,])
#summary(test)

test$安装前最差PH[which(is.na(test$安装前最差PH))]=mean(train$安装前最差PH)#7.278203
test$LAC数值[which(is.na(test$LAC数值))]=mean(train$LAC数值)#6.916143

table(train$event)
table(test$event)
#cox model
items=colnames(train[,-c(1:2)])
s<-paste0(items,collapse = "+")
FML=as.formula(paste0("Surv(Time,event)~",s))
res.cox=coxph(FML,data = train)
summary(res.cox)
# Risk Score = exp(
#   -0.990 × Worst_PH_preECMO +
#     0.0342 × Worst_LAC_preECMO + 
#     0.2056 × Cardiac_arrest_preECMO +
#     -0.6144 × Myocarditis +
#     0.5519 × Bleeding_complication +
#     0.4031 × Renal_complications +
#     0.5489 × Cardiac_complications
# )

#---training set----#
concordance(res.cox,newdata = train)
risk_score=predict(res.cox,newdata = train,type="risk")
roc(train$event,risk_score)
# Data: risk_score in 389 controls (train$event 0) < 224 cases (train$event 1).
# Area under the curve: 0.7434[0.704-0.783]
aa<-roc(train$event,risk_score, plot=TRUE, print.thres=TRUE, ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")

sp.obj <- ci.sp(aa, sensitivities=seq(0, 1, .01), boot.n=100)
plot(sp.obj, type="shape", col="gray60")
#cutoff=1.103 (0.594,0.799)

test_pred=predict(res.cox,newdata = test,type="risk")
roc(test$event,test_pred)
# Data: test_pred in 259 controls (test$event 0) < 149 cases (test$event 1).
# Area under the curve: 0.7299[0.680-0.780]
bb<-roc(test$event,test_pred, plot=TRUE,  ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")

mod.auc=roc(c(train$event,test$event),c(risk_score,test_pred),ci=T)


#timeROC_train
train$lp=predict(res.cox,newdata = train,
                 type = "lp")
time_roc=timeROC(
  T=train$Time,
  delta = train$event,
  marker = train$lp,
  cause = 1,
  weighting = "marginal",
  times = c(3,7,14,28),
  ROC=T,
  iid = T
)
day3<-paste0("3days AUC[95%CI]:",
             sprintf("%.3f",time_roc$AUC[1]),"[",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[1,1]/100),", ",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[1,2]/100),"]")
day7<-paste0("7days AUC[95%CI]:",
             sprintf("%.3f",time_roc$AUC[2]),"[",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[2,1]/100),", ",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[2,2]/100),"]")
day14<-paste0("14days AUC[95%CI]:",
             sprintf("%.3f",time_roc$AUC[3]),"[",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[3,1]/100),", ",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[3,2]/100),"]")
day28<-paste0("28days AUC[95%CI]:",
              sprintf("%.3f",time_roc$AUC[4]),"[",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[4,1]/100),", ",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[4,2]/100),"]")


plot(title="Training Set",time_roc,col = "DodgerBlue",
     time = 3,lty=1,lwd=2,family="serif")
plot(time_roc,col = "LightSeaGreen",add = T,
     time = 7,lty=1,lwd=2,family="serif")
plot(time_roc,col = "DarkOrange",add = T,
     time = 14,lty=1,lwd=2,family="serif")
plot(time_roc,col = "#ED0009",add = T,
     time = 28,lty=1,lwd=2,family="serif")
legend("bottomright",c(day3,day7,day14,day28),
       col=c("DodgerBlue","LightSeaGreen","DarkOrange","#ED0009"),
       lty=1,lwd=2)
#dev.off()

#youden index
train$risk_score=risk_score
roc.obj=roc(train$event,risk_score)
coords(roc.obj,"best",ret="youden")#1.392937


#timeROC_test
test$lp=predict(res.cox,newdata = test,
                 type = "lp")
test_roc=timeROC(
  T=test$Time,
  delta = test$event,
  marker = test$lp,
  cause = 1,
  weighting = "marginal",
  times = c(3,7,14,28),
  ROC=T,
  iid = T
)
day3<-paste0("3day AUC[95%CI]:",
             sprintf("%.3f",test_roc$AUC[1]),"[",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[1,1]/100),", ",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[1,2]/100),"]")
day7<-paste0("7day AUC[95%CI]:",
             sprintf("%.3f",test_roc$AUC[2]),"[",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[2,1]/100),", ",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[2,2]/100),"]")
day14<-paste0("14day AUC[95%CI]:",
              sprintf("%.3f",test_roc$AUC[3]),"[",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[3,1]/100),", ",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[3,2]/100),"]")
day28<-paste0("28days AUC[95%CI]:",
              sprintf("%.3f",test_roc$AUC[4]),"[",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[4,1]/100),", ",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[4,2]/100),"]")

plot(title="Test Set",test_roc,col = "DodgerBlue",
     time = 3,lty=1,lwd=2,family="serif")
plot(test_roc,col = "LightSeaGreen",add = T,
     time = 7,lty=1,lwd=2,family="serif")
plot(test_roc,col = "DarkOrange",add = T,
     time = 14,lty=1,lwd=2,family="serif")
plot(test_roc,col = "#ED0009",add = T,
     time = 28,lty=1,lwd=2,family="serif")
legend("bottomright",c(day3,day7,day14,day28),
       col=c("DodgerBlue","LightSeaGreen","DarkOrange","#ED0009"),
       lty=1,lwd=2)

########train plot ###########
#确定最佳阶段点分高低风险
train$risk_score=risk_score
res.cut=surv_cutpoint(train,time = "Time",
                      event = "event",
                      variables = "risk_score")
plot(res.cut)
cut_point=as.numeric(res.cut$cutpoint[1]) 
train$group="Low"
train$group[which(train$risk_score>cut_point)]="High"
#boxplot death vs alive
colnames(train)
rs.df=train[,c(1,10)]
rs.df$group="death"
rs.df$group[which(rs.df$event==0)]="survival"
ggplot(rs.df,aes(x=group,y=risk_score,fill=group))+
  geom_boxplot(aes(color=group),
               alpha=0.1)+
  geom_jitter(aes(color=group))+
  scale_y_continuous(limits = c(0,10))+
  stat_compare_means(label = "p.format",
                     method = "wilcox.test",
                     label.y = 8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")
  

#KM plot
red="#BC3C28"
blue="#003399"
fit=survfit(Surv(Time,event)~group,data = train)
hist(train$Time)
ggsurvplot(fit,
           pval = T,
           risk.table = T,
           risk.table.col="strata",
           palette = c(red,blue),
           conf.int = T,
           conf.int.alpha=0.1,
           surv.median.line = "hv",
           #xlim=c(0,80),
           break.x.by=20,
           legend.labs=c("High","Low"))
#HR
train$group=as.factor(train$group)
train$group=relevel(train$group,ref = "Low")
summary(coxph(Surv(Time,event)~group,data = train))
#4.0999 [2.955, 5.689] p<2e-16 ***



########test survival plot ###########
rs_test=predict(res.cox,newdata = test,type="risk")
test$risk_score=rs_test
test$group="Low"
test$group[which(test$risk_score>cut_point)]="High"
#boxplot death vs alive
colnames(test)
rs.df=test[,c(1,11)]
rs.df$group="death"
rs.df$group[which(rs.df$event==0)]="survival"
ggplot(rs.df,aes(x=group,y=risk_score,fill=group))+
  geom_boxplot(aes(color=group),
               alpha=0.1)+
  geom_jitter(aes(color=group))+
  scale_y_continuous(limits = c(0,10))+
  stat_compare_means(label = "p.format",
                     method = "wilcox.test",
                     label.y = 8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")

#KM-plot
fit=survfit(Surv(Time,event)~group,data = test)
hist(test$Time)
ggsurvplot(fit,
           pval = T,
           risk.table = T,
           risk.table.col="strata",
           palette = c(red,blue),
           conf.int = T,
           conf.int.alpha=0.1,
           surv.median.line = "hv",
           #xlim=c(0,80),
           break.x.by=20,
           legend.labs=c("High","Low"))
#HR
test$group=as.factor(test$group)
test$group=relevel(test$group,ref = "Low")
summary(coxph(Surv(Time,event)~group,data = test))
#2.7922 [1.91, 4.081] p=1.14e-07 ***

###################
# nomogram plot  #
###################
library(regplot)
library(survival)
summary(train)
train$event=as.numeric(train$event)
colnames(train)[3:9]=c("Worst_PH_preECMO","Worst_LAC_preECMO",
"Cardiac_arrest_preECMO","Myocarditis","Bleeding_complication",
"Renal_complications","Cardiac_complications")


ddist=datadist(train[,1:9])
options(datadist='ddist')
FML <- Surv(Time, event) ~ Worst_PH_preECMO + Worst_LAC_preECMO + 
  Cardiac_arrest_preECMO + Myocarditis + Bleeding_complication + 
  Renal_complications + Cardiac_complications

f_cph=cph(as.formula(FML),data = train,surv = T,time.inc = 3,x = TRUE, y = TRUE)
print(f_cph)


new_pt <- data.frame(
  Time = 2.318056,
  event = 1,
  Worst_PH_preECMO = 7.050,
  Worst_LAC_preECMO = 11.4,
  Cardiac_arrest_preECMO = 0,
  Myocarditis = 0,
  Bleeding_complication = 0,
  Renal_complications = 1,
  Cardiac_complications = 0
)
# 以下是cutoff point的患者数据，对应total point=276
# 0
# 22.802083
# 7.310000
# 10.700000
# 0
# 0
# 0

pn=regplot(f_cph,observation=new_pt, failtime=c(7,14,28),title="Survival Nomogram", 
           prfail=F, points=T,showP = F)

#save(pn,file="./Results/3.nomogramLabel_beforeECMO.RData")
#
load("./Results/3.nomogramLabel_beforeECMO.RData")
