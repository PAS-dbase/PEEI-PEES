setwd("E:/Rworkspace/ECMO_revise")
library(readxl)
library(plyr)
library(rrtable)
library(autoReg)
library(dplyr)
library(data.table)
library(broom)
library(rstatix)
library(forestmodel)#install
library(tidyverse)#install
library(skimr)#install
library(do)#install
library(verification)#install
#patients=read_xlsx("./Data/人口学信息.xlsx",sheet = 1)

##################
# summary table1 #
##################
df=read_xlsx("./Data/人口学信息.xlsx",sheet="表1一般临床资料")
colnames(df)
#ft=gaze(`辅助最终结果`~.,data=df,digits=3)%>%myfit()
summary(df)
#factor
df_f=as.data.frame(lapply(df[,c(1,3,6:8,10:16,21,
                                25:33,36:42,44:48,
                                51:60,62:66,
                                68:72,74:77)],
                          as.factor))
#numeric
df_n=as.data.frame(lapply(df[,c(4,5,17:20)],as.numeric))
df_use=cbind(df_n,df_f)
t1=df_use%>%
  group_by(`辅助最终结果`)%>%
  summarise(across(everything(),~list(summary(.))))
t1=t(t1)
write.csv(t1,"./Results/t1.csv")
#P value
#df_use$ID=rownames(df_use)
table(df_use$支持年份...35)
df_use$支持年份...35=gsub("[^0-9]","",df_use$支持年份...35)
#number compare
df_ln=melt(df_use[,1:7],
           id.var="辅助最终结果",
           measure.var=colnames(df_use)[1:6],
           variable.names="variable",
           value.name = "value")
pn.res=df_ln%>%
  filter(!is.na(value) & !is.na(`辅助最终结果`))%>%
  group_by(variable)%>%
  wilcox_test(data=.,value~`辅助最终结果`)%>%
  add_significance("p")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
#factor comparation
colnames(df_use)
df_lf=melt(df_use[,c(7:64)],
           id.var="辅助最终结果",
           measure.var=colnames(df_use)[8:64],
           variable.names="variable",
           value.name = "value")

summary(df_lf)
#df_lf$value=as.factor(df_lf$value)
pf.res=df_lf%>%
  filter(!is.na(value) & !is.na(`辅助最终结果`))%>%
  group_by(variable)%>%
  kruskal_test(data=.,value~`辅助最终结果`)%>%
  add_significance("p")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
  # summarise(
  #   p.value=chisq.test(table(value,`辅助最终结果`))$p.value
  #   #.group='drop'
  # )
  
p.res=rbind(pn.res[,c(1,7:11)],pf.res[,c(1,4,6,8:10)])
write.csv(p.res,"./Results/pValue_table1.csv",row.names = F)

##################
# summary table2 #
##################
df=read_xlsx("./Data/人口学信息.xlsx",sheet="表2")
colnames(df)
summary(df)
#factor
df_f=as.data.frame(lapply(df[,c(1,9:13,22:24)],
                          as.factor))
#numeric
df_n=as.data.frame(lapply(df[,-c(1,9:13,22:24)],
                          as.numeric))
df_use=cbind(df_n,df_f)
t2=df_use%>%
  group_by(`辅助最终结果`)%>%
  summarise(across(everything(),~list(summary(.))))
t2=t(t2)
write.csv(t2,"./Results/t2.csv")


#P value
#number compare
df_ln=melt(df_use[,1:16],
           id.var="辅助最终结果",
           measure.var=colnames(df_use)[1:15],
           variable.names="variable",
           value.name = "value")
pn.res=df_ln%>%
  filter(!is.na(value) & !is.na(`辅助最终结果`))%>%
  group_by(variable)%>%
  wilcox_test(data=.,value~`辅助最终结果`)%>%
  add_significance("p")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
#factor comparation
colnames(df_use)
df_lf=melt(df_use[,c(16:24)],
           id.var="辅助最终结果",
           measure.var=colnames(df_use)[17:24],
           variable.names="variable",
           value.name = "value")

pf.res=df_lf%>%
  filter(!is.na(value) & !is.na(`辅助最终结果`))%>%
  group_by(variable)%>%
  kruskal_test(data=.,value~`辅助最终结果`)%>%
  add_significance("p")%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

p.res=rbind(pn.res[,c(1,7:11)],pf.res[,c(1,4,6,8:10)])
write.csv(p.res,"./Results/pValue_table2.csv",row.names = F)
