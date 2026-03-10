setwd("E:/Rworkspace/ECMO_revise")
library(dplyr)
library(timeROC)
library(pROC)
library(rms)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(readxl)
library(impute)
library(mice)
library(stringr)
library(ggalluvial)
library(survival)
library(riskRegression)
library(tibble)
library(purrr)
#PEES模型（Prognosis ECMO Evaluation after Support）用于上机24小时后的再次评估
# ==================== part 1 data deal ====================
# 所有差异指标 单cox 
bm_df=read.csv("./Results/2.unicox_data_withoutFill_0.1.csv")
colnames(bm_df)
# ecmo前数据 
data=read.csv("./Results/2.final_Data.csv")
colnames(data)
data$Time=as.numeric(data$Time)

colnames(data)
data=cbind(data,bm_df[,c(10:16,18,20,32,33,37,38)])
summary(data)

data[,17:22]=as.data.frame(lapply(data[,17:22],as.factor))

#模型1 marker
ecmo_pre_vars<-colnames(data)[3:9]
#ecmo后 后补marker
ecmo_post_candidates<-colnames(data)[10:22]

#分训练测试集6：4
death=data[which(data$event==1),]
alive=data[which(data$event==0),]
set.seed(699)
death_random=sample(1:373,224)
alive_random=sample(1:648,389)

#training set
train=rbind(death[death_random,],alive[alive_random,])
#fill NA with knn in the train set
data_imputed <- impute.knn(as.matrix(train[, 1:9]))$data
train<-cbind(as.data.frame(data_imputed),train[,10:22])   
#test set
test=rbind(death[-death_random,],alive[-alive_random,])
#summary(test)

test$安装前最差PH[which(is.na(test$安装前最差PH))]=mean(train$安装前最差PH)
test$LAC数值[which(is.na(test$LAC数值))]=mean(train$LAC数值)

# ==================== part 2 base model ====================
# 基准模型-模型1：仅包含ECMO前变量
base_formula <- as.formula(paste("Surv(Time, event) ~", 
                                 paste(ecmo_pre_vars, collapse = " + ")))

base_model <- coxph(base_formula, data = train,x = TRUE, y = TRUE)
summary(base_model)
saveRDS(base_model, file = "./Results/PEEI_model.rds")
#评估 Base model 的区分度与误差，时间依赖AUC和Brier分数


# 计算基准模型性能
get_model_performance <- function(model, train_data, test_data) {
  # 训练集性能
  c_train <- concordance(model, newdata = train_data)$concordance
  train_risk <- predict(model, newdata = train_data, type = "risk")
  train_auc <- roc(train_data$event, train_risk)$auc
  
  # 测试集性能
  c_test=concordance(model, newdata = test_data)$concordance
  test_risk <- predict(model, newdata = test_data, type = "risk")
  test_auc <- roc(test_data$event, test_risk)$auc
  
  # 汇总
  return(list(
    train_cindex = c_train,
    train_auc = as.numeric(train_auc),
    test_cindex = c_test,
    test_auc = as.numeric(test_auc)
  ))
}

cat("=== 基准模型性能 (仅ECMO前指标) ===\n")
base_perf <- get_model_performance(base_model, train, test)
print(base_perf)

# ==================== part 3：PEES missing value dealing ====================
ecmo_post_candidates#所有ecmo后指标
vars_to_impute=ecmo_post_candidates
# 缺失值多重插补数据集
create_mice_data <- function(data, vars_to_impute) {
  # 选择用于插补的变量（包括结局变量和重要预测变量）
  imp_vars <- c("Time", "event", ecmo_pre_vars, vars_to_impute)
  
  # 创建插补数据集
  imp_data <- data[imp_vars]
  
  # 设置插补方法
  method <- make.method(imp_data)
  #method[imp_data] <- "pmm"
  method["休克分24H"] <- "polyreg"  # 多分类使用polyreg
  method[c("MCS24H","IABP使用","CRRT","经皮穿刺","穿刺和切开")] <- "logreg"  # 二分类使用logreg
  # 进行多重插补（5个数据集）
  mice_imp <- mice(imp_data, m = 5, maxit = 10, method = method, seed = 123)
  
  return(mice_imp)
}

mice_train <- create_mice_data(train, ecmo_post_candidates)
# 1.检查插补迭代过程（收敛诊断）
#plot(mice_train)
# 2. 数值统计比较 
#Mean_Diff：应该小于原始数据标准差的10%
# 提取第一个插补数据集
imp_data_1 <- complete(mice_train, 1)
comparison_results <- data.frame()

for (var in ecmo_post_candidates[1:7]) {
  # 原始数据（去除缺失值）
  original <- na.omit(train[[var]])
  imputed <- imp_data_1[[var]]
  
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
# 24H指标明显更好-----选择24H指标
# Variable Original_Mean Imputed_Mean Original_SD Imputed_SD Mean_Diff
# 1   ECMO4h时PH         7.367        7.370       0.127      0.128     0.003
# 2  ECMO24h时PH         7.397        7.391       0.213      0.228     0.005
# 3  ECMO4h时LAC         5.318        5.037       5.158      4.909     0.281
# 4 ECMO24h时LAC         3.344        3.443       5.446      5.273     0.099
# 5      E4H流量        81.152       81.346      45.594     45.655     0.194
# 6     E12H流量        82.021       82.121      46.132     46.103     0.100
# 7     E24H流量        80.643       80.699      47.257     47.136     0.056
# N_Observed N_Imputed relative_diff
# 1        371       242   0.023622047
# 2        466       147   0.023474178
# 3        363       250   0.054478480
# 4        466       147   0.018178480
# 5        585        28   0.004254946
# 6        584        29   0.002167693
# 7        582        31   0.001185010

mice_test <- create_mice_data(test, ecmo_post_candidates)
test_imp_data <- complete(mice_test, 1)
colnames(test)
colnames(test_imp_data)
test_imp=cbind(test[,c(1:9,11,13,16)],test_imp_data[,c(17:22)])
#fill with mean value
test_imp$ECMO24h时PH[is.na(test_imp$ECMO24h时PH)]=7.391#平均值补齐
test_imp$ECMO24h时LAC[is.na(test_imp$ECMO24h时LAC)]=3.443#平均值补齐
test_imp$E24H流量[is.na(test_imp$E24H流量)]=80.699#平均值补齐
#write.csv(test_imp,"./Results/4.fmodelTestSet_meanFill.csv",row.names = F)
#write.csv(test_imp_data[,c(1:9,13)],"./Results/4.fmodelTestSet.csv",row.names = F)
#write.csv(test_imp_data,"./Results/4.allTestSet.csv",row.names = F)
# ==================== part 4: marker selection ====================
# 逐个纳入最重要的指标
enhanced_models <- list()
colnames(imp_data_1)
# imp_data_1$delta_ph=imp_data_1$ECMO24h时PH-imp_data_1$安装前最差PH
# imp_data_1$delta_lac=imp_data_1$ECMO24h时LAC-imp_data_1$LAC数值
#
# test_imp_data$delta_ph=test_imp_data$ECMO24h时PH-test_imp_data$安装前最差PH
# test_imp_data$delta_lac=test_imp_data$ECMO24h时LAC-test_imp_data$LAC数值
# 从临床最重要的指标开始逐个尝试
recommended_24h_vars <- c("ECMO24h时LAC", "ECMO24h时PH", "E24H流量", 
                          "IABP使用","CRRT","经皮穿刺", "穿刺和切开",
                          "MCS24H", "休克分24H")
summary(imp_data_1)
model_performance <- data.frame()
for (i in seq_along(recommended_24h_vars)) {
  current_var <- recommended_24h_vars[i]
  cat("\n", str_glue("尝试纳入: {current_var}"), "\n")
  
  # 构建新公式
  new_formula <- as.formula(paste("Surv(Time, event) ~", 
                                  paste(c(ecmo_pre_vars, current_var), collapse = " + ")))
  
  # 拟合新模型（添加错误处理）
  new_model <- tryCatch({
    coxph(new_formula, data = imp_data_1)
  }, error = function(e) {
    cat("模型拟合错误:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(new_model)) {
    cat("跳过变量", current_var, "- 模型拟合失败\n")
    next
  }
  # 评估性能（训练集和测试集）
  new_perf <- get_model_performance(new_model, imp_data_1, test_imp)#_data
  
  # 似然比检验
  lrtest <- anova(base_model, new_model)
  lr_pvalue <- lrtest[2, "Pr(>|Chi|)"]
  
  # 安全地获取新变量的系数和显著性
  var_summary <- summary(new_model)$coefficients[8, ]
  var_coef <- var_summary["coef"]
  var_pvalue <- var_summary["Pr(>|z|)"]
  var_hr <- exp(var_coef)
  
  # 存储模型
  enhanced_models[[current_var]] <- new_model
  
  # 记录性能比较
  model_performance <- rbind(model_performance, data.frame(
    Added_Variable = current_var,
    Train_AUC = round(new_perf$train_auc, 4),
    Test_AUC = round(new_perf$test_auc, 4),
    Train_AUC_Improvement = round(new_perf$train_auc - base_perf$train_auc, 4),
    Test_AUC_Improvement = round(new_perf$test_auc - base_perf$test_auc, 4),
    Var_Coefficient = round(var_coef, 4),
    Var_HR = round(var_hr, 4),
    Var_Pvalue = round(var_pvalue, 4),
    LRT_Pvalue = round(lr_pvalue, 4),
    Significant = ifelse(var_pvalue < 0.05, "Yes", "No"),
    Status = "Included"
  ))
  
  # 打印结果
  cat("训练集AUC: ", round(base_perf$train_auc, 4), " -> ", round(new_perf$train_auc, 4), 
      " (提升: +", round(new_perf$train_auc - base_perf$train_auc, 4), ")\n", sep = "")
  cat("测试集AUC: ", round(base_perf$test_auc, 4), " -> ", round(new_perf$test_auc, 4), 
      " (提升: +", round(new_perf$test_auc - base_perf$test_auc, 4), ")\n", sep = "")
  cat("变量系数:", round(var_coef, 4), "HR:", round(var_hr, 4), 
      "p-value:", round(var_pvalue, 4), "\n")
  cat("似然比检验 p-value:", round(lr_pvalue, 4), "\n")
}

#mean impute
# 尝试纳入: ECMO24h时LAC 
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# 训练集AUC: 0.7434 -> 0.7596 (提升: +0.0162)
# 测试集AUC: 0.7299 -> 0.7425 (提升: +0.0126)
# 变量系数: 0.0376 HR: 1.0383 p-value: 0 
# 似然比检验 p-value: 0 
# 
# 尝试纳入: ECMO24h时PH 
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# 训练集AUC: 0.7434 -> 0.7511 (提升: +0.0077)
# 测试集AUC: 0.7299 -> 0.7364 (提升: +0.0065)
# 变量系数: -0.517 HR: 0.5963 p-value: 0.0038 
# 似然比检验 p-value: 0.0155 
# 
# 尝试纳入: E24H流量 
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# 训练集AUC: 0.7434 -> 0.7437 (提升: +2e-04)
# 测试集AUC: 0.7299 -> 0.7353 (提升: +0.0054)
# 变量系数: 0.0017 HR: 1.0017 p-value: 0.2027 
# 似然比检验 p-value: 0.2244 

# 尝试纳入: MCS24H 
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# 训练集AUC: 0.7434 -> 0.7457 (提升: +0.0022)
# 测试集AUC: 0.7299 -> 0.7315 (提升: +0.0016)
# 变量系数: -0.446 HR: 0.6402 p-value: 0.2518 
# 似然比检验 p-value: 0.2203 
# 
# 尝试纳入: 休克分24H 
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# Setting levels: control = 0, case = 1
# Setting direction: controls < cases
# 训练集AUC: 0.7434 -> 0.7618 (提升: +0.0184)
# 测试集AUC: 0.7299 -> 0.7333 (提升: +0.0034)
# 变量系数: 0.3111 HR: 1.3649 p-value: 0.255 
# 似然比检验 p-value: 0 

#逐步纳入图形展示
plot_data <- model_performance %>%
  dplyr::mutate(
    Variable = factor(Added_Variable, levels = rev(Added_Variable)),
    Significance = case_when(
      Var_Pvalue < 0.001 ~ "***",
      Var_Pvalue < 0.01 ~ "**",
      Var_Pvalue < 0.05 ~ "*",
      Var_Pvalue < 0.1 ~ ".",
      TRUE ~ "NS"
    ),
    Label = sprintf("+%.3f (p=%s)", Test_AUC_Improvement, Significant),
    Color = ifelse(Var_Pvalue < 0.05, "#2E8B57", "#8B0000")  # 显著绿色，不显著红色
  )

p1 <- ggplot(plot_data, aes(x = Test_AUC_Improvement, y = reorder(Variable, Test_AUC_Improvement))) +
  geom_segment(aes(xend = 0, yend = Variable), color = "gray50") +
  geom_point(aes(color = Var_Pvalue < 0.05, size = abs(Test_AUC_Improvement)), 
             alpha = 0.8) +
  geom_text(aes(label = sprintf("Δ=+%.3f\np=%.4f", Test_AUC_Improvement, Var_Pvalue),
                x = Test_AUC_Improvement + 0.003),
            size = 3, hjust = 0, vjust = 0.5) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#27AE60"),
                     labels = c("not sig", "sig")) +
  labs(
    title = "after-ECMO marker contribution",
    x = "AUC improve in the test set",
    y = "",
    color = "significance",
    size = "improvement"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  ) +
  guides(size = "none")  # 隐藏size图例

print(p1)

included_vars <- model_performance$Added_Variable[model_performance$Status == "Included"]
# ==================== part 5: final model ====================
significant_vars <- c("ECMO24h时LAC", "ECMO24h时PH")

# 方案A：只纳入ECMO24h时LAC（效果最好的单个变量）
final_formula_single <- as.formula(paste("Surv(Time, event) ~", 
                                         paste(c(ecmo_pre_vars, "ECMO24h时LAC"), collapse = " + ")))

# 方案B：纳入所有显著变量
final_formula_combined <- as.formula(paste("Surv(Time, event) ~", 
                                           paste(c(ecmo_pre_vars, significant_vars), collapse = " + ")))

# 拟合两个最终模型
final_model_single <- coxph(final_formula_single, data = imp_data_1)
final_model_combined <- coxph(final_formula_combined, data = imp_data_1)
summary(final_model_single)
summary(final_model_combined)
# 比较性能
final_perf_single <- get_model_performance(final_model_single, imp_data_1, test_imp)#_data
final_perf_combined <- get_model_performance(final_model_combined, imp_data_1, test_imp)#_data

cat("\n=== 最终模型性能比较 ===\n")
cat("1. 仅ECMO24h时LAC:\n")
cat("   训练集AUC:", round(final_perf_single$train_auc, 4), 
    "测试集AUC:", round(final_perf_single$test_auc, 4), "\n")

cat("2. ECMO24h时LAC + ECMO24h时PH:\n")
cat("   训练集AUC:", round(final_perf_combined$train_auc, 4), 
    "测试集AUC:", round(final_perf_combined$test_auc, 4), "\n")
# # 选择最佳最终模型
# if (final_perf_combined$test_auc > final_perf_single$test_auc) {
#   final_model <- final_model_combined
#   final_perf <- final_perf_combined
#   cat("\n*** 选择组合模型: ECMO24h时LAC + ECMO24h时PH ***\n")
# } else {
#   final_model <- final_model_single
#   final_perf <- final_perf_combined
#   cat("\n*** 选择单个模型: 仅ECMO24h时LAC ***\n")
# }
#--------模型比较图---------#
#training set
train_risk_base <- predict(base_model, newdata = imp_data_1, type = "risk")
train_risk_single <- predict(final_model_single, newdata = imp_data_1, type = "risk")
train_risk_combined <- predict(final_model_combined, newdata = imp_data_1, type = "risk")

train_roc_base <- roc(imp_data_1$event, train_risk_base,plot=TRUE,  ci=TRUE,
                      print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
train_roc_single <- roc(imp_data_1$event, train_risk_single,plot=TRUE,  ci=TRUE,
                        print.auc=TRUE,legacy.axes = TRUE,col = "#C783CB",add=T)
train_roc_combined <- roc(imp_data_1$event, train_risk_combined,plot=TRUE,  ci=TRUE,
                          print.auc=TRUE,legacy.axes = TRUE,col = "#F58F5B",add=T)

legend("bottomright", 
       legend = c(paste("base model (AUC =", round(train_roc_base$auc, 3), ")"),
                  paste("base model+24hLAC (AUC =", round(train_roc_single$auc, 3), ")"),
                  paste("base model+24hLAC+24hPH (AUC =", round(train_roc_combined$auc, 3), ")")),
       col = c("#1687A7", "#C783CB", "#F58F5B"), lwd = 2, cex = 0.8)

#test set
test_risk_base <- predict(base_model, newdata = test_imp, type = "risk")#_data
test_risk_single <- predict(final_model_single, newdata = test_imp, type = "risk")
test_risk_combined <- predict(final_model_combined, newdata = test_imp, type = "risk")

test_roc_base <- roc(test_imp_data$event, test_risk_base,plot=TRUE,  ci=TRUE,
                      print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
test_roc_single <- roc(test_imp_data$event, test_risk_single,plot=TRUE,  ci=TRUE,
                        print.auc=TRUE,legacy.axes = TRUE,col = "#C783CB",add=T)
test_roc_combined <- roc(test_imp_data$event, test_risk_combined,plot=TRUE,  ci=TRUE,
                          print.auc=TRUE,legacy.axes = TRUE,col = "#F58F5B",add=T)


legend("bottomright", 
       legend = c(paste("base model (AUC =", round(test_roc_base$auc, 3), ")"),
                  paste("base model+24hLAC (AUC =", round(test_roc_single$auc, 3), ")"),
                  paste("base model+24hLAC+24hPH (AUC =", round(test_roc_combined$auc, 3), ")")),
       col = c("#1687A7", "#C783CB", "#F58F5B"), lwd = 2, cex = 0.8)

# 敏感性分析：检查pH在不同亚组中的表现
# 按乳酸水平分组，看pH是否有额外价值
# length(imp_data_1$ECMO24h时LAC[which(imp_data_1$ECMO24h时LAC>4)])
# imp_data_1$lac_group <- ifelse(imp_data_1$ECMO24h时LAC > 4, "高乳酸", "低乳酸")
# test_imp_data$lac_group <- ifelse(test_imp_data$ECMO24h时LAC > 4, "高乳酸", "低乳酸")
# for (group in unique(imp_data_1$lac_group)) {
#   subgroup <- imp_data_1[imp_data_1$lac_group == group, ]
#   if (nrow(subgroup) > 50) {
#     model_sub <- coxph(Surv(Time, event) ~ ECMO24h时PH, data = subgroup)
#     p_value <- summary(model_sub)$coefficients["ECMO24h时PH", "Pr(>|z|)"]
#     cat("在", group, "组中，pH的p-value:", round(p_value, 4), "\n")
#   }
# }
# model_sub1=coxph(final_formula_combined, data = imp_data_1[which(imp_data_1$lac_group=="低乳酸"),])
# summary(model_sub1)

#final model
# 拟合最终模型
final_model <- coxph(final_formula_single, data = imp_data_1)
saveRDS(final_model, file = "./Results/PEES_model.rds")
# 最终模型性能
final_perf <- get_model_performance(final_model, imp_data_1, test_imp)#_data

cat("\n=== 最终组合模型性能 ===\n")
cat("训练集AUC:", round(final_perf$train_auc, 4), 
    "(基准:", round(base_perf$train_auc, 4), 
    ", 提升: +", round(final_perf$train_auc - base_perf$train_auc, 4), ")\n")
cat("测试集AUC:", round(final_perf$test_auc, 4), 
    "(基准:", round(base_perf$test_auc, 4), 
    ", 提升: +", round(final_perf$test_auc - base_perf$test_auc, 4), ")\n")

# 最终模型摘要
summary(final_model)
# Risk Score = exp(
#   -0.958 × Worst_PH_preECMO +
#     0.029 × Worst_LAC_preECMO + 
#     0.207 × Cardiac_arrest_preECMO +
#     -0.631 × Myocarditis +
#     0.571 × Bleeding_complication +
#     0.388 × Renal_complications +
#     0.526 × Cardiac_complications+
#     0.038×ECMO24LAC
# )

cat("=== ECMO24h时LAC临床意义说明 ===\n\n")
lac_hr <- 1.0383
calculate_risk_changes <- function(lac_reduction) {
  risk_reduction <- (1 - (1 / (lac_hr ^ lac_reduction))) * 100
  return(round(risk_reduction, 1))
}

cat("🔬 具体临床意义：\n")
cat("• 乳酸每下降 1 mmol/L → 死亡风险降低", calculate_risk_changes(1), "%\n")
# ==================== plot ====================
imp_data_1$base_risk_score <- predict(base_model, newdata = imp_data_1, type = "risk")
test_imp$base_risk_score <- predict(base_model, newdata = test_imp, type = "risk")#_data


# 按基准模型分层
imp_data_1$base_risk_group <- ifelse(imp_data_1$base_risk_score > 1.103, 
                                     "base_high", "base_low")
test_imp$base_risk_group <- ifelse(test_imp$base_risk_score > 1.103, 
                                     "base_high", "base_low")



########train plot ###########
risk_score_train=predict(final_model,newdata = imp_data_1,type="risk")
imp_data_1$risk_score=risk_score_train
#time ROC
time_roc=timeROC(
  T=imp_data_1$Time,
  delta = imp_data_1$event,
  marker = imp_data_1$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(7,14,28),
  ROC=T,
  iid = T
)
day7<-paste0("7days AUC[95%CI]:",
             sprintf("%.3f",time_roc$AUC[1]),"[",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[1,1]/100),", ",
             sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[1,2]/100),"]")
day14<-paste0("14days AUC[95%CI]:",
              sprintf("%.3f",time_roc$AUC[2]),"[",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[2,1]/100),", ",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[2,2]/100),"]")
day28<-paste0("28days AUC[95%CI]:",
              sprintf("%.3f",time_roc$AUC[3]),"[",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[3,1]/100),", ",
              sprintf("%.3f",confint(time_roc,level=0.95)$CI_AUC[3,2]/100),"]")
plot(title="Training Set",time_roc,col = "DodgerBlue",
     time = 7,lty=1,lwd=2,family="serif")
plot(time_roc,col = "LightSeaGreen",add = T,
     time = 14,lty=1,lwd=2,family="serif")
plot(time_roc,col = "DarkOrange",add = T,
     time = 28,lty=1,lwd=2,family="serif")
legend("bottomright",c(day7,day14,day28),
       col=c("DodgerBlue","LightSeaGreen","DarkOrange"),
       lty=1,lwd=2)


#确定最佳阶段点分高低风险
res.cut=surv_cutpoint(imp_data_1,time = "Time",
                      event = "event",
                      variables = "risk_score")
plot(res.cut)# 1.574694
cut_point=as.numeric(res.cut$cutpoint[1]) 
imp_data_1$final_risk_group="final_low"
imp_data_1$final_risk_group[which(imp_data_1$risk_score>cut_point)]="final_high"
table(imp_data_1$final_risk_group)
#boxplot death vs alive
colnames(imp_data_1)
rs.df=imp_data_1[,c(2,19)]
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
fit=survfit(Surv(Time,event)~final_risk_group,data = imp_data_1)
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
imp_data_1$final_risk_group=as.factor(imp_data_1$final_risk_group)
imp_data_1$final_risk_group=relevel(imp_data_1$final_risk_group,ref = "final_low")
summary(coxph(Surv(Time,event)~final_risk_group,data = imp_data_1))
#3.628  [ 2.786，4.725]   p<2e-16 ***

########test plot ###########
#time ROC
risk_score_test=predict(final_model,newdata = test_imp,type="risk")#_data
test_imp$risk_score=risk_score_test
test_roc=timeROC(
  T=test_imp$Time,
  delta = test_imp$event,
  marker = test_imp$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(7,14,28),
  ROC=T,
  iid = T
)
day7<-paste0("7days AUC[95%CI]:",
             sprintf("%.3f",test_roc$AUC[1]),"[",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[1,1]/100),", ",
             sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[1,2]/100),"]")
day14<-paste0("14days AUC[95%CI]:",
              sprintf("%.3f",test_roc$AUC[2]),"[",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[2,1]/100),", ",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[2,2]/100),"]")
day28<-paste0("28days AUC[95%CI]:",
              sprintf("%.3f",test_roc$AUC[3]),"[",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[3,1]/100),", ",
              sprintf("%.3f",confint(test_roc,level=0.95)$CI_AUC[3,2]/100),"]")
plot(title="Test Set",test_roc,col = "DodgerBlue",
     time = 7,lty=1,lwd=2,family="serif")
plot(test_roc,col = "LightSeaGreen",add = T,
     time = 14,lty=1,lwd=2,family="serif")
plot(test_roc,col = "DarkOrange",add = T,
     time = 28,lty=1,lwd=2,family="serif")
legend("bottomright",c(day7,day14,day28),
       col=c("DodgerBlue","LightSeaGreen","DarkOrange"),
       lty=1,lwd=2)

#确定最佳阶段点分高低风险
test_imp$final_risk_group="final_low"
test_imp$final_risk_group[which(test_imp$risk_score>cut_point)]="final_high"
table(test_imp$final_risk_group)
#boxplot death vs alive
colnames(test_imp)
rs.df=test_imp[,c(1,17)]
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
fit=survfit(Surv(Time,event)~final_risk_group,data = test_imp)
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
test_imp$final_risk_group=as.factor(test_imp$final_risk_group)
test_imp$final_risk_group=relevel(test_imp$final_risk_group,ref = "final_low")
summary(coxph(Surv(Time,event)~final_risk_group,data = test_imp))
#3.023  [2.218,4.188]   p=2.88e-11 ***

# ==========评估LAC是否可以作为可靠的单一指标========
# 首先确定24h乳酸的cutoff（用于着色）
roc_24hlac <- roc(imp_data_1$event, imp_data_1$ECMO24h时LAC)
lac_cutoff <- coords(roc_24hlac, "best", ret = "threshold")$threshold

cat("24h乳酸最佳cutoff:", round(lac_cutoff, 2), "mmol/L\n")
#24h乳酸最佳cutoff: 2.55 mmol/L


#按照指南大于2是高乳酸血症
lac_cutoff=2
# 创建重分类数据，并添加乳酸信息
reclass_data <- rbind(imp_data_1[,c(13,19:22)],test_imp[,c(11,15:18)]) %>%
  mutate(
    # 重分类分组 - 使用准确的列值
    reclassification = case_when(
      base_risk_group == "base_low" & final_risk_group == "final_low" ~ "Low",
      base_risk_group == "base_high" & final_risk_group == "final_high" ~ "High", 
      base_risk_group == "base_low" & final_risk_group == "final_high" ~ "Worse",
      base_risk_group == "base_high" & final_risk_group == "final_low" ~ "Recover"
    ),
    # 乳酸风险分组
    lac_risk = ifelse(ECMO24h时LAC > lac_cutoff, "Lac_high", "Lac_low")
  )
reclass_data$base_risk_group=factor(reclass_data$base_risk_group)
reclass_data$final_risk_group=factor(reclass_data$final_risk_group)
reclass_data$lac_risk=factor(reclass_data$lac_risk)
reclass_data$event=c(imp_data_1$event,test_imp$event)
write.csv(reclass_data,"./Results/4.reclass_data.csv",row.names = F)

class_colors <- c(
  "High"    = "#D55E00",
  "Low"     = "#0072B2",
  "Recover" = "#01A089",
  "Worse"   = "#E54C33"
)

df <- ddply(
  reclass_data,
  .(base_risk_group, final_risk_group,reclassification, lac_risk),
  summarise,
  Freq = length(reclassification)
)
# 总人数
N = sum(df$Freq)

# 计算百分比
df$percent <- df$Freq / N * 100
df$label <- sprintf("%.1f%%", df$percent)
write.csv(df,"./Results/4.classification.csv",row.names = F)

palluvial=ggplot(df,
                 aes(axis1 = base_risk_group,
                     axis2 = lac_risk,
                     axis3 = final_risk_group,
                     y = Freq,
                     fill = reclassification)) +
  geom_alluvium(alpha = 0.85, width = 1/12) +
  geom_stratum(fill = "grey92", color = "grey40", width = 1/5) +
  geom_label(stat = "stratum",
             aes(label = after_stat(stratum)),
             size = 4.3, label.size = 0,
             fill = "white", alpha = 0.9) +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(limits = c("Base risk", "24h LAC", "Final classification"),
                   expand = c(.05, .05)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold")
  ) +
  ggtitle("ECMO Reclassification Pathway\n(Base → 24h LAC → Final)")

df$Freq
print(palluvial)
147/(235+101+84+147)#0.2592593
# reclass_table <- table(
#   reclassification = reclass_data$reclassification,
#   lac_risk = reclass_data$lac_risk
# )
# # 转换为数据框并计算比例
# reclass_summary <- as.data.frame(reclass_table)
# reclass_summary <- reclass_summary %>%
#   group_by(reclassification) %>%
#   mutate(
#     total = sum(Freq),
#     proportion = Freq / total
#   )
# 
# 计算最终模型风险与乳酸风险的一致性
consistency_table <- table(
  最终模型风险 = reclass_data$final_risk_group,
  乳酸风险 = ifelse(reclass_data$ECMO24h时LAC > lac_cutoff, "高风险", "低风险")
)

cat("=== 最终模型风险 vs 乳酸风险一致性分析 ===\n")
print(consistency_table)

# 计算Kappa一致性系数
library(vcd)
kappa_result <- Kappa(consistency_table)
# value    ASE      z  Pr(>|z|)
# Unweighted -0.286 0.0282 -10.14 3.653e-24
# Weighted   -0.286 0.0282 -10.14 3.653e-24
# 计算一致率
agreement_rate <- sum(diag(consistency_table)) / sum(consistency_table) * 100
cat("总体一致率:", round(agreement_rate, 1), "%\n")
#总体一致率: 34.9 %

# 计算结局
ct=table(
  结局 = reclass_data$event,
  乳酸风险 = ifelse(reclass_data$ECMO24h时LAC > lac_cutoff, "高风险", "低风险")
)
kappa_ct <- Kappa(ct)
# value     ASE     z  Pr(>|z|)
# Unweighted 0.2083 0.02986 6.975 3.052e-12
# Weighted   0.2083 0.02986 6.975 3.052e-12
# ==========评估LAC清除率是否可以作为可靠的单一指标======
colnames(imp_data_1)
colnames(test_imp_data)
# ΔLAC = 安装前 - 24H
imp_data_1$delta_lac <- imp_data_1$LAC数值 - imp_data_1$ECMO24h时LAC
test_imp_data$delta_lac <- test_imp_data$LAC数值 - test_imp_data$ECMO24h时LAC

# 清除率 = (安装前 - 24H) / 安装前
imp_data_1$lac_clearance <- (imp_data_1$LAC数值 - imp_data_1$ECMO24h时LAC) / imp_data_1$LAC数值 * 100
test_imp_data$lac_clearance <- (test_imp_data$LAC数值 - test_imp_data$ECMO24h时LAC) / test_imp_data$LAC数值 * 100


#清除率 cutoff
df_down <- imp_data_1[imp_data_1$delta_lac > 0, ]  # 乳酸下降者

roc_clr_down <- roc(df_down$event, df_down$delta_lac)
roc(df_down$event, df_down$delta_lac, plot=TRUE, print.thres=TRUE, ci=TRUE,
    print.auc=TRUE,legacy.axes = TRUE,col = "#D55E00")
coords_clr_down <- coords(roc_clr_down, "best",
                          ret = c("threshold", "sensitivity", "specificity"),
                          transpose = FALSE)

coords_clr_down

best_cutoff <- coords_clr_down["threshold"]$threshold
# 三类分类
# imp_data_1$lac_resp3="Clear_low"
# imp_data_1$lac_resp3[which(imp_data_1$delta_lac<0)]="Lactate_up"
# imp_data_1$lac_resp3[which(imp_data_1$delta_lac>= best_cutoff)]="Clear_high"

imp_data_1$lac_resp3<- with(imp_data_1, ifelse(
  ECMO24h时LAC > LAC数值,
  "Lactate_up",     # 乳酸上升 → 反应最差
  ifelse(delta_lac >= best_cutoff,
         "Clear_high",    # ΔLAC 大于 cutoff → 清除好
         "Clear_low"      # ΔLAC 小于 cutoff → 清除不足
  )
))

imp_data_1$lac_resp3 <- factor(
  imp_data_1$lac_resp3,
  levels = c("Lactate_up", "Clear_low", "Clear_high")
)
table(imp_data_1$lac_resp3)

#绘图
df_sankey <- ddply(
  imp_data_1,
  .(base_risk_group, lac_resp3, final_risk_group),
  summarise,
  Freq = length(lac_resp3)
)
lac_colors <- c(
  "Lactate_up" = "#01A089",   # 红：恶化
  "Clear_low"  = "#4BBBD3",   # 灰：轻度改善
  "Clear_high" = "#E54C33"    # 绿：明显改善
)

pSankey=ggplot(df_sankey,
               aes(axis1 = base_risk_group,
                   axis2 = lac_resp3,
                   axis3 = final_risk_group,
                   y = Freq,
                   fill = lac_resp3)) +   
  
  geom_alluvium(alpha = 0.85, width = 1/12) +
  geom_stratum(fill = "grey95", color = "grey40", width = 1/6) +
  geom_label(stat = "stratum",
             aes(label = after_stat(stratum)),
             size = 4.2, label.size = 0,
             fill = "white", alpha = 0.9) +
  
  scale_fill_manual(values = lac_colors) +
  scale_x_discrete(
    limits = c("Base risk", "ΔLactate response", "Final risk"),
    expand = c(.05, .05)
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  
  ggtitle("ECMO Response Pathway Visualization\n(Base → ΔLactate response → Final risk)")

print(pSankey)



###################
# nomogram plot  #
###################
library(regplot)
library(survival)
summary(imp_data_1)
traindf=imp_data_1[,c(1:9,13)]
traindf$event=as.numeric(traindf$event)


traindf=read.csv("./Results/4.fmodelTrainingSet.csv")
final_model<- readRDS("./Results/EPES_model.rds")
risk_score_train=predict(final_model,newdata = traindf,type="risk")
traindf$risk_score=risk_score_train

colnames(traindf)[3:10]=c("Worst_PH_preECMO","Worst_LAC_preECMO",
                       "Cardiac_arrest_preECMO","Myocarditis","Bleeding_complication",
                       "Renal_complications","Cardiac_complications","LAC_24h")

ddist=datadist(traindf[,1:10])
options(datadist='ddist')
FML <- Surv(Time, event) ~ Worst_PH_preECMO + Worst_LAC_preECMO + 
  Cardiac_arrest_preECMO + Myocarditis + Bleeding_complication + 
  Renal_complications + Cardiac_complications+LAC_24h

f_cph=cph(as.formula(FML),data = traindf,surv = T,time.inc = 3,x = TRUE, y = TRUE)
print(f_cph)

# med=Quantile(f_cph)
# surv=Survival(f_cph)
cut_pt <- data.frame(
  Time = 1.898611,
  event = 1,
  Worst_PH_preECMO = 7.12,
  Worst_LAC_preECMO = 1.5,
  Cardiac_arrest_preECMO = 0,
  Myocarditis = 0,
  Bleeding_complication = 0,
  Renal_complications = 1,
  Cardiac_complications = 0,
  LAC_24h=5.3
)#total point=130


pn=regplot(f_cph, observation=cut_pt, failtime=c(7,14,28),
           title="Survival Nomogram", prfail=F, clickable=F, points=TRUE)#cut patient

new_pt <- data.frame(
  Time = 2.318056,
  event = 1,
  Worst_PH_preECMO = 7.050,
  Worst_LAC_preECMO = 11.4,
  Cardiac_arrest_preECMO = 0,
  Myocarditis = 0,
  Bleeding_complication = 0,
  Renal_complications = 1,
  Cardiac_complications = 0,
  LAC_24h=5.9
)

pn=regplot(f_cph, observation=new_pt, failtime=c(7,14,28),
           title="Survival Nomogram", prfail=F, clickable=F, points=TRUE)

save(pn,file="./Results/4.nomogramLabel_afterECMO.RData")





