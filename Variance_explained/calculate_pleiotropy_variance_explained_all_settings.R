# =============================================================================
# 计算多效性SNP对表型解释量的百分比（所有四种设置）
# BothPos, BothNeg, AposBneg, AnegBpos
# =============================================================================
rm(list=ls())

# 模拟参数
N_Ind_Train <- 20000
N_Eff_SNP_PhenoA <- 100
N_Eff_SNP_PhenoB <- 100

# 不同设置对应的随机种子（与原始模拟代码保持一致）
# BothPos: seed=73, BothNeg: seed=52, AposBneg: seed=52, AnegBpos: seed=36
setting_seeds <- list(
  "BothPos" = 73,
  "BothNeg" = 52,
  "AposBneg" = 52,
  "AnegBpos" = 36
)

# =============================================================================
# 定义计算函数（可处理不同的多效性效应方向）
# =============================================================================

calculate_variance_decomposition <- function(setting_name, pleio_effect, n_pleio, atob_val, precomputed_SNP_effects) {
  # 根据设置确定多效性SNP的效应方向
  if (setting_name == "BothPos") {
    pleio_mean_A <- pleio_effect
    pleio_mean_B <- pleio_effect
  } else if (setting_name == "BothNeg") {
    pleio_mean_A <- -pleio_effect
    pleio_mean_B <- -pleio_effect
  } else if (setting_name == "AposBneg") {
    pleio_mean_A <- pleio_effect
    pleio_mean_B <- -pleio_effect
  } else if (setting_name == "AnegBpos") {
    pleio_mean_A <- -pleio_effect
    pleio_mean_B <- pleio_effect
  } else {
    stop("Unknown setting: ", setting_name)
  }
  
  # 使用预先计算的SNP效应大小（确保同一设置下所有参数组合使用相同的SNP效应）
  SNP_Effect_A <- precomputed_SNP_effects$SNP_Effect_A
  SNP_Effect_B <- precomputed_SNP_effects$SNP_Effect_B
  
  # 根据设置使用对应的随机种子（用于生成基因型数据和多效性SNP效应）
  current_seed <- setting_seeds[[setting_name]]
  if (is.null(current_seed)) {
    stop("Unknown setting: ", setting_name)
  }
  set.seed(current_seed)
  
  # 生成基因型数据
  # 注意：为了确保多效性SNP效应的一致性，我们需要使用固定的随机数序列
  # 但由于n_pleio不同，我们需要为每个n_pleio生成不同的多效性SNP效应
  # 因此，我们使用一个基于n_pleio的偏移量来确保可重复性
  Data_full <- matrix(NA, nrow=N_Ind_Train, ncol=(N_Eff_SNP_PhenoA + n_pleio + N_Eff_SNP_PhenoB))
  for (i in 1:ncol(Data_full)) {
    Data_full[,i] <- sample(0:1, N_Ind_Train, replace=TRUE)
  }
  
  # 生成多效性SNP效应（这些会根据n_pleio和pleio_effect变化）
  SNP_Pleiotropy_Effect_A <- rnorm(n_pleio, mean=pleio_mean_A, sd=pleio_effect)
  SNP_Pleiotropy_Effect_B <- rnorm(n_pleio, mean=pleio_mean_B, sd=pleio_effect)
  
  # 计算表型A
  PhenoA_specific <- Data_full[,1:N_Eff_SNP_PhenoA] %*% SNP_Effect_A
  PhenoA_genetic <- PhenoA_specific + 
                    Data_full[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+n_pleio)] %*% SNP_Pleiotropy_Effect_A
  PhenoA_pleio <- Data_full[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+n_pleio)] %*% SNP_Pleiotropy_Effect_A
  var_A_specific <- var(PhenoA_specific)
  var_A_pleio <- var(PhenoA_pleio)
  var_A_genetic <- var(PhenoA_genetic)
  PhenoA <- PhenoA_genetic + rnorm(length(PhenoA_genetic), mean=0, sd=sd(PhenoA_genetic)*2)
  var_A_total <- var(PhenoA)
  var_A_env <- var_A_total - var_A_genetic
  
  # 计算表型B
  PhenoB_causal <- PhenoA_genetic * atob_val
  PhenoB_pleio <- Data_full[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+n_pleio)] %*% SNP_Pleiotropy_Effect_B
  PhenoB_specific <- Data_full[,(N_Eff_SNP_PhenoA+n_pleio+1):(N_Eff_SNP_PhenoA+n_pleio+N_Eff_SNP_PhenoB)] %*% SNP_Effect_B
  PhenoB_genetic <- PhenoB_causal + PhenoB_pleio + PhenoB_specific
  var_B_genetic <- var(PhenoB_genetic)
  PhenoB <- PhenoB_genetic + rnorm(length(PhenoB_genetic), mean=0, sd=sd(PhenoB_genetic)*2)
  var_B_total <- var(PhenoB)
  var_B_env <- var_B_total - var_B_genetic
  
  # 计算百分比（表型A - 由于各成分独立，可以直接计算）
  pct_A_specific_genetic <- var_A_specific / var_A_genetic * 100
  pct_A_specific_total <- var_A_specific / var_A_total * 100
  pct_A_pleio_genetic <- var_A_pleio / var_A_genetic * 100
  pct_A_pleio_total <- var_A_pleio / var_A_total * 100
  pct_A_env_total <- var_A_env / var_A_total * 100
  
  # 计算百分比（表型B - 使用协方差方法确保总和=100%）
  # 各成分对总遗传方差的贡献 = cov(总遗传, 各成分) / var(总遗传)
  # 这样可以确保总和等于100%
  cov_B_causal_genetic <- cov(PhenoB_genetic, PhenoB_causal)
  cov_B_pleio_genetic <- cov(PhenoB_genetic, PhenoB_pleio)
  cov_B_specific_genetic <- cov(PhenoB_genetic, PhenoB_specific)
  
  # 占遗传方差百分比（考虑协方差，总和=100%）
  pct_B_causal_genetic <- cov_B_causal_genetic / var_B_genetic * 100
  pct_B_pleio_genetic <- cov_B_pleio_genetic / var_B_genetic * 100
  pct_B_specific_genetic <- cov_B_specific_genetic / var_B_genetic * 100
  
  # 占总表型方差百分比（使用协方差方法）
  cov_B_causal_total <- cov(PhenoB, PhenoB_causal)
  cov_B_pleio_total <- cov(PhenoB, PhenoB_pleio)
  cov_B_specific_total <- cov(PhenoB, PhenoB_specific)
  
  pct_B_causal_total <- cov_B_causal_total / var_B_total * 100
  pct_B_pleio_total <- cov_B_pleio_total / var_B_total * 100
  pct_B_specific_total <- cov_B_specific_total / var_B_total * 100
  pct_B_env_total <- var_B_env / var_B_total * 100
  
  return(data.frame(
    Setting = setting_name,
    Pleio_Effect = pleio_effect,
    N_Pleio_SNP = n_pleio,
    AtoB = atob_val,
    # 表型A
    PhenoA_Specific_Pct_Genetic = round(pct_A_specific_genetic, 2),
    PhenoA_Specific_Pct_Total = round(pct_A_specific_total, 2),
    PhenoA_Pleio_Pct_Genetic = round(pct_A_pleio_genetic, 2),
    PhenoA_Pleio_Pct_Total = round(pct_A_pleio_total, 2),
    PhenoA_Env_Pct_Total = round(pct_A_env_total, 2),
    # 表型B
    PhenoB_Causal_Pct_Genetic = round(pct_B_causal_genetic, 2),
    PhenoB_Causal_Pct_Total = round(pct_B_causal_total, 2),
    PhenoB_Pleio_Pct_Genetic = round(pct_B_pleio_genetic, 2),
    PhenoB_Pleio_Pct_Total = round(pct_B_pleio_total, 2),
    PhenoB_Specific_Pct_Genetic = round(pct_B_specific_genetic, 2),
    PhenoB_Specific_Pct_Total = round(pct_B_specific_total, 2),
    PhenoB_Env_Pct_Total = round(pct_B_env_total, 2)
  ))
}

# =============================================================================
# 预先为每个设置生成SNP效应大小（确保同一设置下所有参数组合使用相同的SNP效应）
# =============================================================================

cat(rep("=", 80), "\n")
cat("预先生成SNP效应大小（确保同一设置下的一致性）\n")
cat(rep("=", 80), "\n\n")

# 为每个设置预先生成SNP效应大小
SNP_effects_by_setting <- list()
for (setting in c("BothPos", "BothNeg", "AposBneg", "AnegBpos")) {
  current_seed <- setting_seeds[[setting]]
  set.seed(current_seed)
  
  # 生成SNP效应大小（与原始模拟代码保持一致）
  SNP_Effect_A <- rnorm(N_Eff_SNP_PhenoA, mean=2, sd=2)
  SNP_Effect_B <- rnorm(N_Eff_SNP_PhenoB, mean=2, sd=2)
  
  SNP_effects_by_setting[[setting]] <- list(
    SNP_Effect_A = SNP_Effect_A,
    SNP_Effect_B = SNP_Effect_B
  )
  cat("已为设置", setting, "生成SNP效应大小\n")
}

# =============================================================================
# 不同参数组合的计算
# =============================================================================

cat(rep("=", 80), "\n")
cat("不同设置下的表型方差分解分析\n")
cat(rep("=", 80), "\n\n")

settings <- c("BothPos", "BothNeg", "AposBneg", "AnegBpos")
pleio_effects <- c(0.5, 1, 1.5, 2)
pleio_snp_counts <- c(10, 20, 30, 40, 50)
AtoB_values <- c(0, 0.1, 0.2, 0.3)

results_df <- data.frame()

for (setting in settings) {
  cat("正在计算设置:", setting, "\n")
  for (pleio_eff in pleio_effects) {
    for (n_pleio in pleio_snp_counts) {
      for (atob_val in AtoB_values) {
        result <- calculate_variance_decomposition(setting, pleio_eff, n_pleio, atob_val, 
                                                   SNP_effects_by_setting[[setting]])
        results_df <- rbind(results_df, result)
      }
    }
  }
}

# 保存结果
result_file <- "pleiotropy_variance_explained_all_settings.csv"
write.csv(results_df, result_file, row.names = FALSE)
cat("\n结果已保存到:", result_file, "\n")
cat("总行数:", nrow(results_df), "\n\n")

# =============================================================================
# 显示汇总统计
# =============================================================================

cat("汇总统计（Pleiotropic Effect = 1, AtoB = 0）：\n\n")

for (setting in settings) {
  subset_data <- subset(results_df, Setting == setting & Pleio_Effect == 1 & AtoB == 0)
  cat("设置:", setting, "\n")
  cat("  表型A - 占遗传方差百分比：\n")
  cat("    特异性SNP - 均值:", round(mean(subset_data$PhenoA_Specific_Pct_Genetic), 2), "%\n")
  cat("    多效性SNP - 均值:", round(mean(subset_data$PhenoA_Pleio_Pct_Genetic), 2), "%\n")
  cat("  表型B - 占遗传方差百分比：\n")
  cat("    A->B因果效应 - 均值:", round(mean(subset_data$PhenoB_Causal_Pct_Genetic), 2), "%\n")
  cat("    多效性SNP - 均值:", round(mean(subset_data$PhenoB_Pleio_Pct_Genetic), 2), "%\n")
  cat("    特异性SNP - 均值:", round(mean(subset_data$PhenoB_Specific_Pct_Genetic), 2), "%\n\n")
}

cat("汇总统计（Pleiotropic Effect = 1, N_Pleio_SNP = 30, 不同AtoB值）：\n\n")

for (setting in settings) {
  subset_data <- subset(results_df, Setting == setting & Pleio_Effect == 1 & N_Pleio_SNP == 30)
  cat("设置:", setting, "\n")
  for (atob in AtoB_values) {
    row <- subset_data[subset_data$AtoB == atob, ]
    if (nrow(row) > 0) {
      cat(sprintf("  AtoB = %.1f: A->B因果效应 = %.2f%%, 多效性SNP = %.2f%%, 特异性SNP = %.2f%%\n",
                  atob, row$PhenoB_Causal_Pct_Genetic, row$PhenoB_Pleio_Pct_Genetic, 
                  row$PhenoB_Specific_Pct_Genetic))
    }
  }
  cat("\n")
}

