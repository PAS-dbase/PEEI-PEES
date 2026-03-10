setwd("E:/Rworkspace/ECMO_revise")
#########selected marker without fill######
df=read.csv("./Results/1.cleaning_data_withoutFill.csv")
##############
# univar cox #
##############
library(dplyr)
library(survival)
library(survminer)
library(plyr)
library(MASS)
library(rms)
library(splines)
library(glmnet)
library(riskRegression)
library(timeROC)
colnames(df)
# #PH、lac等变化情况
# df$delta_PH=df$ECMO24h时PH-df$安装前最差PH
# #乳酸清除量
# df$delta_LAC=df$LAC数值-df$ECMO24h时LAC
# #血管活性药物评分减少值
# df$delta_VIS=df$VIS数值-df$E24小时VIS 
# write.csv(df,"./Results/2.cleaning_data_withoutFill_delta.csv",row.names = F)
#delete imbalance variables
# 检查每个候选二分类/因子变量的分布
check_tab <- function(varname){
  tab <- table(df[[varname]], df$event)
  print(tab)
}
check_tab("经人工血管")#只有两例数据经人工血管
################################
# remove imbalance variables  #
###############################
#######---自动筛选：剔除样本数太少的变量---#####
factor_vars <- colnames(df)[23:69]
summary(df$年龄)
# 筛选条件：
# - 每个水平至少10人
# - 每个水平至少有5个事件 & 5个非事件
valid_vars <- c()

for (v in factor_vars) {
  tab <- table(df[[v]], df$event)
  # 跳过空表或全缺失
  if (nrow(tab) == 0 || ncol(tab) < 2) next
  sizes <- rowSums(tab)
  cond1 <- all(sizes >= 10)
  cond2 <- all(tab[, "0"] >= 5 & tab[, "1"] >= 5)  # 事件列名可能是0/1，确认一下
  if (cond1 && cond2) {
    valid_vars <- c(valid_vars, v)
  } else {
    message("变量 ", v, " 被剔除（不平衡）")
    print(tab)
  }
}
#######结果#########
# 变量 心律失常 被剔除（不平衡）
# 
#    0   1
# 0 633 371
# 1  15   2
# 变量 外伤 被剔除（不平衡）
# 
#    0   1
# 0 644 361
# 1   4  12
# 变量 代谢紊乱 被剔除（不平衡）
# 
#    0   1
# 0 641 371
# 1   7   2
# 变量 长期口服抗凝药物史 被剔除（不平衡）
# 
#     0   1
# 0 507 305
# 1   2   4
# 变量 其他地点 被剔除（不平衡）
# 
#    0   1
# 0 643 373
# 1   5   0
# 变量 经人工血管 被剔除（不平衡）
# 
#    0   1
# 0 648 371
# 1   0   2

valid_vars
setdiff(factor_vars,valid_vars)#
# [1] "ECMO前漂浮导管检查" "心律失常"           "外伤"              
# [4] "代谢紊乱"           "长期口服抗凝药物史" "其他地点"          
# [7] "经人工血管"
df_del_imbalance=df[,!colnames(df)%in%setdiff(factor_vars,valid_vars)]
write.csv(df_del_imbalance,"./Results/2.cleaning_data_deleteImbalance.csv",row.names = F)

df=read.csv("./Results/2.cleaning_data_deleteImbalance.csv")
summary(df)
surv=Surv(time=df$Time,event =as.numeric(df$event) )
df$surv=surv
colnames(df)
#factor
df[,c(6,23:62)]=as.data.frame(lapply(df[,c(6,23:62)],as.factor))

#unicox
UniCox <- function(x) {
  FML <- as.formula(paste0("surv~", x))
  Cox <- coxph(FML, data = df)
  Sum <- summary(Cox)
  
  coef_df <- as.data.frame(Sum$coefficients)
  confint_df <- as.data.frame(Sum$conf.int)
  
  if (is.factor(df[[x]])) {
    var_levels <- levels(df[[x]])
    ref_level <- var_levels[1]  
    
    if (nrow(coef_df) > 1) {
      current_rows <- rownames(coef_df)
      compare_levels <- sub(paste0("^", x), "", current_rows)
      
      var_names <- paste0(x, ": ", compare_levels, " vs ", ref_level)
    } else {
      var_names <- x
    }
  } else {
    if (nrow(coef_df) > 1) {
      var_names <- paste0(x, "_", rownames(coef_df))
    } else {
      var_names <- x
    }
  }
  
  HR <- round(coef_df[, "exp(coef)"], 3)
  Pval <- round(coef_df[, "Pr(>|z|)"], 3)
  
  if (nrow(confint_df) > 0) {
    CI_lower <- round(confint_df[, "lower .95"], 3)
    CI_upper <- round(confint_df[, "upper .95"], 3)
    CI <- paste0(CI_lower, "-", CI_upper)
  } else {
    CI <- rep(NA, nrow(coef_df))
  }
  
  Unicox <- data.frame(
    "variables" = var_names,
    "HR" = HR,
    "95%CI" = CI,
    "Pvalue" = Pval,
    stringsAsFactors = FALSE
  )
  
  return(Unicox)
}

varNames=colnames(df)[6:65]
UniVar=list()
for (nm in varNames) {
  temp=tryCatch(
    {UniCox(nm)},
    error=function(e){message('Error @',nm);return(NULL)},
    finally = {message('next')}
  )
  if (!is.null(temp)) {
    UniVar[[nm]] <- temp
  }
}
UniVar=ldply(UniVar,data.frame)
UniVar=UniVar[,2:5]
write.csv(UniVar,"./Results/2.unicox_withoutFill.csv",row.names = F)
writexl::write_xlsx(UniVar,"./Results/2.unicox_withoutFill.xlsx")
# Variables with p < 0.1 in univariate analysis 
# were entered into the multivariable Cox regression.
cox_items=UniVar$varables[which(UniVar$Pvalue<0.1)]

cox_data=df[,which(colnames(df)%in%cox_items)]
cox_data=cbind(df[,1:2],cox_data)
write.csv(cox_data,"./Results/2.unicox_data_withoutFill_0.1.csv",row.names = F)

################################
# before ecmo marker selection #
################################
bm_df=read.csv("./Results/2.unicox_data_withoutFill_0.1.csv")
colnames(bm_df)
bm_df_sel=bm_df[,c(1,2,4:9,17:34)]
str(bm_df_sel)
bm_df_sel$event=as.numeric(bm_df_sel$event)
df_sel_na=na.omit(bm_df_sel)


items=colnames(df_sel_na)[-c(1:2)]
s<-paste0(items,collapse = "+")
FML=as.formula(paste0("Surv(Time,event)~",s))
####--coxph--####
res.cox=coxph(FML,data = df_sel_na)
res.step=stepAIC(res.cox,direction = "both")
mucox=summary(res.step)
mu_items=names(which(mucox$coefficients[,5]<0.05))
mu_items
####--lasso cox---#####
X<-model.matrix(~. -1,data=df_sel_na[,-c(1,2)])
y<-Surv(df_sel_na$Time,df_sel_na$event)
set.seed(123)
cvfit <- cv.glmnet(X,y, family="cox", alpha=1)  # alpha=1=LASSO
fit_lasso <- glmnet(X, y, family = "cox", alpha = 1, lambda = cvfit$lambda.1se)
# 入模变量
sel_lasso <- rownames(coef(fit_lasso))[as.numeric(coef(fit_lasso)) != 0]
sel_lasso

#######################
# method comparison   #
#######################
# Stepwise Cox 风险评分
lp_step <- predict(res.step, type="lp")

# time
tau <- 28   

# timeROC for stepwise
roc_step <- timeROC(T=df_sel_na$Time, delta=df_sel_na$event,
                    marker=lp_step, cause=1, times=tau, iid=TRUE)

lp_lasso <- as.numeric(predict(fit_lasso, newx=X, type="link"))
# timeROC for lasso
roc_lasso <- timeROC(T=df_sel_na$Time, delta=df_sel_na$event,
                     marker=lp_lasso, cause=1, times=tau, iid=TRUE)

cat("Stepwise Cox - 28天AUC:", roc_step$AUC, "\n")
#Stepwise Cox - 28天AUC: NA 0.7622503 
cat("LASSO Cox   - 28天AUC:", roc_lasso$AUC, "\n")
#LASSO Cox   - 28天AUC: NA 0.7181571 
sc_auc <- Score(
  object    = list(step = lp_step, lasso = lp_lasso),  # 两个“风险分数”
  formula   = Surv(Time, event) ~ 1,
  data      = df_sel_na,
  times     = 28,
  metrics   = "AUC",    # 只算 AUC，避免调用 C.survival
  riskScore = TRUE,
  summary   = "ipa")
sc_auc$AUC

# C-index（Harrell's C）
c_step  <- concordance(Surv(Time, event) ~ lp_step,  data = df_sel_na)$concordance
c_lasso <- concordance(Surv(Time, event) ~ lp_lasso, data = df_sel_na)$concordance
c_step; c_lasso

#coxph better than lasso

###############
# multi cox   #
###############
summary(res.step)

#forest plot
library(forestmodel)
pdf("./Figs/2.forestplot.pdf",width = 12,height = 8,family = "GB1")
forest_model(res.step,#cox模型
             theme = theme_forest(),
             factor_separate_line=TRUE
)
dev.off()

data_markers=bm_df[,which(colnames(bm_df)%in%mu_items)]
data_markers=cbind(bm_df[,1:2],data_markers)
write.csv(data_markers,"./Results/2.data_markerSelection.csv",row.names = F)


##########################
# RF importance ranking #
##########################
library(randomForestSRC)
library(ggplot2)
cox_data=na.omit(bm_df_sel)

items=colnames(cox_data)[-c(1:2)]
s<-paste0(items,collapse = "+")
FML=as.formula(paste0("Surv(Time,event)~",s))
rsf_fit=rfsrc(FML,data=cox_data,ntree=1000,
              importance = T)
var_imp=rsf_fit$importance
#plot(get.tree(rsf_fit,3))
var_imp=as.data.frame(var_imp)
var_imp$variables=rownames(var_imp)
var_imp=var_imp[order(-var_imp$var_imp),]
var_imp$var_imp=var_imp$var_imp*(100/sum(var_imp$var_imp))
pdf("./Figs/2.importance_cleaning.pdf",width = 8,height = 6,family = "GB1")
ggplot(data=var_imp,mapping = aes(x=reorder(variables,var_imp),y=var_imp,
                                  fill=variables,group=factor(1)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=round(var_imp,3),color="black"))+
  coord_flip()+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid = element_blank())
dev.off()
imp_items=var_imp$variables[1:10]
#intersect(mu_items,imp_items)
union(mu_items,imp_items)
#final marker selection
final_data=data_markers[,-9]
final_data$`LAC数值`=bm_df_sel$LAC数值
final_data=final_data[,c(1:3,9,4:8)]
write.csv(final_data,"./Results/2.final_Data.csv",row.names = F)

