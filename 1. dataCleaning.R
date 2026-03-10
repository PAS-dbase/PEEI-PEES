setwd("E:/Rworkspace/ECMO_revise")
library(readxl)
library(lubridate)
library(mice)
library(do)
library(autoReg)
library(rrtable)
########################
# Data standardization #
########################
data=read_xlsx("./Data/清洗模型指标+PECSOS评分.xlsx",sheet=2)
table(data$辅助最终结果)
#Time=Discharge time - Ecmo start time
df=data.frame(event=data$辅助最终结果)
di=ymd_hm(data$出院时间)-ymd_hm(data$ECMO开始时间)
df$Time=as.double(di,units="days")
hosday=ymd_hm(data$出院时间)-ymd_hm(data$入院时间)
df$`住院时间`=as.double(hosday,units="days")
ecmotime=ymd_hm(data$ECMO撤机时间)-ymd_hm(data$ECMO开始时间)
df$`ECMO使用时间`=as.double(ecmotime,units="days")
icuday=ymd_hm(data$出ICU时间)-ymd_hm(data$入ICU时间)
df$`ICU停留时间`=as.double(icuday,units="days")
df$`ICU停留时间`[is.na(df$`ICU停留时间`)]=0#没有入ICU的时间为0
colnames(data)
df_use=cbind(df,data[,8:78])
df_use$event[which(df_use$event=="2")]=0 #2=survival 1=dead

#1=no 2=yes
df_use[,c(24:69,71,72,74)]=lapply(df_use[,c(24:69,71,72,74)],function(x) ifelse(x==1,0,x))
df_use[,c(24:69,71,72,74)]=lapply(df_use[,c(24:69,71,72,74)],function(x) ifelse(x==2,1,x))
#factor
df_use[,c(6,24:76)]=as.data.frame(lapply(df_use[,c(6,24:76)],as.factor))
#numeric
table(df_use$ECMO24h时LAC)
df_use$ECMO24h时LAC[which(df_use$ECMO24h时LAC=="＜1")]=0.5
df_use$ECMO24h时LAC[which(df_use$ECMO24h时LAC=="＞15")]=16
df_use$ECMO24h时LAC[which(df_use$ECMO24h时LAC=="＞20")]=21
df_use$ECMO24h时LAC[which(df_use$ECMO24h时LAC=="nd")]=0
df_use$ECMO24h时LAC[which(df_use$ECMO24h时LAC=="ND")]=0
df_use[,7:23]=as.data.frame(lapply(df_use[,7:23],as.numeric))
summary(df_use)
table(df_use$event)
colnames(df_use)[c(20:22,50)]=c("E4H流量","E12H流量","E24H流量","感染并发症")
write.csv(df_use,"./Results/1.standard_data_time.csv",row.names = F)

########################
#  Deal missing data   #
########################
df_use=read.csv("./Results/1.standard_data_time.csv")
# wrong data check and delete
df_del=df_use[-which(df_use$Time<=0),]#外院转运，入院前已经上ecmo
df_del=df_del[-which(df_del$ECMO使用时间<1),] #1021
#delete col which missing value over 50%
summary(df_del)
df_del=df_del[,colSums(is.na(df_del))<511]
#df_del=df_del[,-c(3,5)]
write.csv(df_del,"./Results/1.cleaning_data_withoutFill.csv",row.names = F)

#table1
df_del=read.csv("./Results/1.cleaning_data_withoutFill.csv")
colnames(df_del)
mean_sd(df_del$年龄)
# y      ymin    ymax
# 1 6.718903 0.8297035 12.6081
table(df_del$性别)#gender
table(df_del$event)
df_del[,c(6,23:69)]=as.data.frame(lapply(df_del[,c(6,23:69)],as.factor))
ft <- gaze(event~.,data=df_del,digits = 3) %>%myft()
table2docx(ft,target="./Results/Table1.docx")
# summary
ft <- gaze(event~.,data=df_del,digits = 3,method=4) %>%myft()
table2docx(ft,target="./Results/Table1_iqr.docx")

#Other instead using rf mean values 
colnames(df_del)
imp=mice(df_del[,1:22],
         method = "rf",
         seed=1,
         printFlag = F)  
df_imputed=complete(imp)
df_del[,1:22]=df_imputed
write.csv(df_del,"./Results/1.cleaning_data_withFill.csv",row.names = F)

