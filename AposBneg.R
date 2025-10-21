setwd('/home1/WangXF/MRsimulation/AposBneg/')

rm(list=ls())  # 清理工作环境

# ===================== 文件夹创建函数 =====================
# 创建参数对应的子文件夹
create_parameter_folder <- function(N_Ind_Train, N_Ind_Test, AtoB, pleio_effect, N_Eff_SNP_PhenoAB) {
  folder_name <- paste0("N_Train_", N_Ind_Train, "_",
                        "N_Test_", N_Ind_Test, "_",
                        "AtoB_", AtoB, "_",
                        "Pleio_Effect_", pleio_effect, "_",
                        "Pleio_SNP_", N_Eff_SNP_PhenoAB)
  
  # 创建主文件夹
  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  # 创建测试集子文件夹
  test_folder <- file.path(folder_name, "test_sets")
  if (!dir.exists(test_folder)) {
    dir.create(test_folder, recursive = TRUE)
  }
  
  return(list(main_folder = folder_name, test_folder = test_folder))
}

# ===================== Complete Sensitivity Analysis Simulation =====================

# 主模拟函数：运行PRS和IVW方法的统计效力分析
run_simulation <- function(N_Ind_Train, N_Ind_Test, AtoB, pleio_effect, N_Eff_SNP_PhenoAB, N_Sim=1000, seed=52, folder_path=NULL, test_folder_path=NULL) {
  set.seed(seed)  # 设置随机种子保证可复现
  
  # 基础参数设置
  N_SNP <- 10000  # 总SNP数量
  N_Eff_SNP_PhenoA <- 100  # 影响表型A的SNP数量
  N_Eff_SNP_PhenoB <- 100  # 影响表型B的SNP数量
  N_Threshold <- 10  # PRS阈值数量
  NPerm <- 1000  # permutation次数
  
  # 生成训练集基因型数据（0/1编码）
  Data <- matrix(NA, nrow=N_Ind_Train, ncol=N_SNP)
  for (i in 1:N_SNP) {
    Data[,i] <- sample(0:1, N_Ind_Train, replace=TRUE)  # 随机生成0/1基因型
  }
  
  # 生成SNP效应大小
  SNP_Effect_A <- rnorm(N_Eff_SNP_PhenoA, mean=2, sd=2)  # 影响表型A的SNP效应
  SNP_Pleiotropy_Effect_A <- rnorm(N_Eff_SNP_PhenoAB, mean=pleio_effect, sd=pleio_effect)  # 多效性SNP对表型A的效应
  SNP_Pleiotropy_Effect_B <- rnorm(N_Eff_SNP_PhenoAB, mean=-pleio_effect, sd=pleio_effect)  # 多效性SNP对表型B的效应（相同）
  SNP_Effect_B <- rnorm(N_Eff_SNP_PhenoB, mean=2, sd=2)  # 影响表型B的SNP效应
  
  # 生成训练集表型数据
  PhenoA <- Data[,1:N_Eff_SNP_PhenoA] %*% SNP_Effect_A +  # 表型A：SNP效应 + 多效性效应
    Data[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB)] %*% SNP_Pleiotropy_Effect_A
  PhenoA <- PhenoA + rnorm(length(PhenoA), mean=0, sd=sd(PhenoA)*2)  # 添加环境噪声
  
  PhenoB <- PhenoA*AtoB +  # 表型B：因果效应 + 多效性效应 + SNP效应
    Data[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB)] %*% SNP_Pleiotropy_Effect_B +
    Data[,(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB+N_Eff_SNP_PhenoB)] %*% SNP_Effect_B
  PhenoB <- PhenoB + rnorm(length(PhenoB), mean=0, sd=sd(PhenoB)*2)  # 添加环境噪声
  
  # 训练集GWAS分析
  Covar_A <- cov(Data, PhenoA)  # SNP与表型A的协方差
  Var_SNP <- apply(Data,2,var)  # SNP方差
  Beta_A <- Covar_A/Var_SNP  # SNP对表型A的效应大小
  Cor_A <- cor(Data, PhenoA)  # SNP与表型A的相关性
  T_A <- Cor_A/sqrt(1-Cor_A^2)*sqrt(N_Ind_Train-2)  # t统计量
  P_A <- 1-2*abs(pt(T_A, df=N_Ind_Train-2, lower.tail=F)-0.5)  # P值
  List_A <- P_A<0.05  # 显著SNP标记
  Order_A <- order(P_A)  # P值排序
  
  Covar_B <- cov(Data, PhenoB)  # SNP与表型B的协方差
  Beta_B <- Covar_B/Var_SNP  # SNP对表型B的效应大小
  Cor_B <- cor(Data, PhenoB)  # SNP与表型B的相关性
  T_B <- Cor_B/sqrt(1-Cor_B^2)*sqrt(N_Ind_Train-2)  # t统计量
  P_B <- 1-2*abs(pt(T_B, df=N_Ind_Train-2, lower.tail=F)-0.5)  # P值
  List_B <- P_B<0.05  # 显著SNP标记
  Order_B <- order(P_B)  # P值排序
  
  # PRS SNP选择矩阵（用于构建PRS）
  Keep_A <- matrix(NA, nrow=N_SNP, ncol=N_Threshold)  # 表型A的SNP选择矩阵
  for (i in 1:N_Threshold) {
    Keep_A[,i] <- rep(TRUE,N_SNP)  # 初始全选
    Keep_A[!List_A,i] <- FALSE  # 排除不显著的SNP
    Keep_A[Order_B[1:(i*50)],i] <- FALSE  # 排除与表型B最相关的SNP
  }
  Keep_B <- matrix(NA, nrow=N_SNP, ncol=N_Threshold)  # 表型B的SNP选择矩阵
  for (i in 1:N_Threshold) {
    Keep_B[,i] <- rep(TRUE,N_SNP)  # 初始全选
    Keep_B[!List_B,i] <- FALSE  # 排除不显著的SNP
    Keep_B[Order_A[1:(i*50)],i] <- FALSE  # 排除与表型A最相关的SNP
  }
  
  # 训练集IVW分析（A→B方向，不删除多效性SNP）
  Model_IVW_AtoB_NoFilter_Train <- lm(Beta_B[List_A] ~ -1 + Beta_A[List_A])  # 使用所有显著SNP
  IVW_AtoB_NoFilter_Estimate_Train <- summary(Model_IVW_AtoB_NoFilter_Train)$coef[1,1]  # 因果效应估计
  IVW_AtoB_NoFilter_P_Train <- summary(Model_IVW_AtoB_NoFilter_Train)$coef[1,4]  # P值
  
  # 训练集IVW分析（A→B方向，删除多效性SNP）
  Model_IVW_AtoB_Filter_Train <- lm(Beta_B[List_A & !List_B] ~ -1 + Beta_A[List_A & !List_B])  # 排除同时显著的SNP
  IVW_AtoB_Filter_Estimate_Train <- summary(Model_IVW_AtoB_Filter_Train)$coef[1,1]  # 因果效应估计
  IVW_AtoB_Filter_P_Train <- summary(Model_IVW_AtoB_Filter_Train)$coef[1,4]  # P值
  
  # 训练集IVW分析（B→A方向，不删除多效性SNP）
  Model_IVW_BtoA_NoFilter_Train <- lm(Beta_A[List_B] ~ -1 + Beta_B[List_B])  # 使用所有显著SNP
  IVW_BtoA_NoFilter_Estimate_Train <- summary(Model_IVW_BtoA_NoFilter_Train)$coef[1,1]  # 因果效应估计
  IVW_BtoA_NoFilter_P_Train <- summary(Model_IVW_BtoA_NoFilter_Train)$coef[1,4]  # P值
  
  # 训练集IVW分析（B→A方向，删除多效性SNP）
  Model_IVW_BtoA_Filter_Train <- lm(Beta_A[List_B & !List_A] ~ -1 + Beta_B[List_B & !List_A])  # 排除同时显著的SNP
  IVW_BtoA_Filter_Estimate_Train <- summary(Model_IVW_BtoA_Filter_Train)$coef[1,1]  # 因果效应估计
  IVW_BtoA_Filter_P_Train <- summary(Model_IVW_BtoA_Filter_Train)$coef[1,4]  # P值
  
  # 测试集模拟结果存储
  P_Sim_AtoB_Overall <- rep(0,N_Sim)  # PRS A→B方法P值
  P_Sim_BtoA_Overall <- rep(0,N_Sim)  # PRS B→A方法P值
  P_Sim_IVW_AtoB_NoFilter_Test <- rep(0,N_Sim)  # IVW A→B不删SNP方法P值
  P_Sim_IVW_AtoB_Filter_Test <- rep(0,N_Sim)  # IVW A→B删SNP方法P值
  P_Sim_IVW_BtoA_NoFilter_Test <- rep(0,N_Sim)  # IVW B→A不删SNP方法P值
  P_Sim_IVW_BtoA_Filter_Test <- rep(0,N_Sim)  # IVW B→A删SNP方法P值
  
  # PRS工具变量验证结果存储
  PRS_IV_cor_A <- rep(0,N_Sim)  # PRS A与表型A相关性
  PRS_IV_Pperms_A <- rep(0,N_Sim)  # PRS A permutation p值
  PRS_IV_cor_B <- rep(0,N_Sim)  # PRS B与表型B相关性
  PRS_IV_Pperms_B <- rep(0,N_Sim)  # PRS B permutation p值
  
  # 开始测试集模拟 - 使用并行计算
  sim_results <- foreach(k = 1:N_Sim, .combine = 'rbind', .packages = c('stats'), 
                         .export = c('test_folder_path')) %dopar% {
                           # 生成测试集基因型数据
                           Data_2 <- matrix(NA, nrow=N_Ind_Test, ncol=N_SNP)
                           for (i in 1:N_SNP) {
                             Data_2[,i] <- sample(0:1, N_Ind_Test, replace=TRUE)  # 随机生成0/1基因型
                           }
                           
                           # 生成测试集表型数据
                           PhenoA2 <- Data_2[,1:N_Eff_SNP_PhenoA] %*% SNP_Effect_A +  # 表型A：SNP效应 + 多效性效应
                             Data_2[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB)] %*% SNP_Pleiotropy_Effect_A
                           PhenoA2 <- PhenoA2 + rnorm(length(PhenoA2), mean=0, sd=sd(PhenoA2)*2)  # 添加环境噪声
                           
                           PhenoB2 <- PhenoA2*AtoB +  # 表型B：因果效应 + 多效性效应 + SNP效应
                             Data_2[,(N_Eff_SNP_PhenoA+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB)] %*% SNP_Pleiotropy_Effect_B +
                             Data_2[,(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB+1):(N_Eff_SNP_PhenoA+N_Eff_SNP_PhenoAB+N_Eff_SNP_PhenoB)] %*% SNP_Effect_B
                           PhenoB2 <- PhenoB2 + rnorm(length(PhenoB2), mean=0, sd=sd(PhenoB2)*2)  # 添加环境噪声
                           
                           # 保存当前测试集的基因型和表型数据
                           test_sim_data <- list(
                             Data_Test = Data_2,
                             PhenoA_Test = PhenoA2,
                             PhenoB_Test = PhenoB2,
                             Sim_Index = k,
                             N_Ind_Test = N_Ind_Test,
                             N_SNP = N_SNP,
                             AtoB = AtoB,
                             Pleio_Effect = pleio_effect,
                             N_Eff_SNP_PhenoAB = N_Eff_SNP_PhenoAB
                           )
                           
                           # 为每个测试集生成唯一文件名并保存到测试集子文件夹
                           test_filename <- paste0("test_data_sim_", k, "_N", N_Ind_Test, "_AtoB", AtoB, "_Pleio", pleio_effect, "_SNP", N_Eff_SNP_PhenoAB, ".RData")
                           test_filepath <- file.path(test_folder_path, test_filename)
                           save(test_sim_data, file = test_filepath)
                           
                           # PRS方法分析（A→B方向）
                           PRS_A <- matrix(NA, nrow=N_Ind_Test, ncol=N_Threshold)  # PRS A矩阵
                           for (i in 1:N_Threshold) {
                             PRS_A[,i] <- Data_2[,Keep_A[,i]] %*% Beta_A[Keep_A[,i]]  # 构建PRS A
                           }
                           Cor_AtoB <- cor(PRS_A, PhenoB2)  # PRS A与表型B的相关性
                           Cor_AtoB_Perm <- matrix(NA, nrow=ncol(PRS_A), ncol=NPerm)  # permutation相关性矩阵
                           for (i in 1:NPerm) {
                             Perm <- sample(1:length(PhenoB2),replace=F)  # 随机打乱个体顺序
                             Cor_AtoB_Perm[,i] <- cor(PRS_A[Perm,],PhenoB2)  # 计算permutation相关性
                           }
                           P_Sim_AtoB_Overall_k <- mean(abs(mean(Cor_AtoB)) < abs(colMeans(Cor_AtoB_Perm)))  # permutation p值
                           
                           # PRS方法分析（B→A方向）
                           PRS_B <- matrix(NA, nrow=N_Ind_Test, ncol=N_Threshold)  # PRS B矩阵
                           for (i in 1:N_Threshold) {
                             PRS_B[,i] <- Data_2[,Keep_B[,i]] %*% Beta_B[Keep_B[,i]]  # 构建PRS B
                           }
                           Cor_BtoA <- cor(PRS_B, PhenoA2)  # PRS B与表型A的相关性
                           Cor_BtoA_Perm <- matrix(NA, nrow=ncol(PRS_B), ncol=NPerm)  # permutation相关性矩阵
                           for (i in 1:NPerm) {
                             Perm <- sample(1:length(PhenoA2),replace=F)  # 随机打乱个体顺序
                             Cor_BtoA_Perm[,i] <- cor(PRS_B[Perm,],PhenoA2)  # 计算permutation相关性
                           }
                           P_Sim_BtoA_Overall_k <- mean(abs(mean(Cor_BtoA)) < abs(colMeans(Cor_BtoA_Perm)))  # permutation p值
                           
                           # PRS工具变量验证（A→B方向）
                           PRS_Cor_with_PhenoA <- cor(PRS_A, PhenoA2)  # PRS A与表型A的相关性
                           PRS_IV_cor_A_k <- mean(PRS_Cor_with_PhenoA)  # 平均相关性
                           
                           # PRS工具变量permutation验证（A→B方向）
                           PRS_IV_Perm_Cor_AtoB <- matrix(NA, nrow=ncol(PRS_A), ncol=NPerm)  # permutation相关性矩阵
                           for (perm in 1:NPerm) {
                             Perm <- sample(1:length(PhenoA2),replace=F)  # 随机打乱个体顺序
                             PRS_Cor_perm <- cor(PRS_A[Perm,],PhenoA2)  # 计算permutation相关性
                             PRS_IV_Perm_Cor_AtoB[,perm] <- PRS_Cor_perm
                           }
                           PRS_IV_Pperms_A_k <- mean(abs(mean(PRS_Cor_with_PhenoA)) < abs(colMeans(PRS_IV_Perm_Cor_AtoB)))  # permutation p值
                           
                           # PRS工具变量验证（B→A方向）
                           PRS_Cor_with_PhenoB <- cor(PRS_B, PhenoB2)  # PRS B与表型B的相关性
                           PRS_IV_cor_B_k <- mean(PRS_Cor_with_PhenoB)  # 平均相关性
                           
                           # PRS工具变量permutation验证（B→A方向）
                           PRS_IV_Perm_Cor_BtoA <- matrix(NA, nrow=ncol(PRS_B), ncol=NPerm)  # permutation相关性矩阵
                           for (perm in 1:NPerm) {
                             Perm <- sample(1:length(PhenoB2),replace=F)  # 随机打乱个体顺序
                             PRS_Cor_perm <- cor(PRS_B[Perm,],PhenoB2)  # 计算permutation相关性
                             PRS_IV_Perm_Cor_BtoA[,perm] <- PRS_Cor_perm
                           }
                           PRS_IV_Pperms_B_k <- mean(abs(mean(PRS_Cor_with_PhenoB)) < abs(colMeans(PRS_IV_Perm_Cor_BtoA)))  # permutation p值
                           
                           # 测试集GWAS分析
                           Covar_A2 <- cov(Data_2, PhenoA2)  # SNP与表型A的协方差
                           Covar_B2 <- cov(Data_2, PhenoB2)  # SNP与表型B的协方差
                           Var_SNP2 <- apply(Data_2,2,var)  # SNP方差
                           Beta_A2 <- Covar_A2/Var_SNP2  # SNP对表型A的效应大小
                           Beta_B2 <- Covar_B2/Var_SNP2  # SNP对表型B的效应大小
                           
                           # 测试集IVW分析（A→B方向，不删除多效性SNP）
                           P_Sim_IVW_AtoB_NoFilter_Test_k <- 1  # 默认值
                           tryCatch({
                             Model_IVW_AtoB_NoFilter_Test <- lm(Beta_B2[List_A] ~ -1 + Beta_A2[List_A])  # 使用所有显著SNP
                             P_Sim_IVW_AtoB_NoFilter_Test_k <- summary(Model_IVW_AtoB_NoFilter_Test)$coef[1,4]  # P值
                           }, error = function(e) {
                             P_Sim_IVW_AtoB_NoFilter_Test_k <- 1  # 出错时设为1
                           })
                           
                           # 测试集IVW分析（A→B方向，删除多效性SNP）
                           P_Sim_IVW_AtoB_Filter_Test_k <- 1  # 默认值
                           tryCatch({
                             Model_IVW_AtoB_Filter_Test <- lm(Beta_B2[List_A & !List_B] ~ -1 + Beta_A2[List_A & !List_B])  # 排除同时显著的SNP
                             P_Sim_IVW_AtoB_Filter_Test_k <- summary(Model_IVW_AtoB_Filter_Test)$coef[1,4]  # P值
                           }, error = function(e) {
                             P_Sim_IVW_AtoB_Filter_Test_k <- 1  # 出错时设为1
                           })
                           
                           # 测试集IVW分析（B→A方向，不删除多效性SNP）
                           P_Sim_IVW_BtoA_NoFilter_Test_k <- 1  # 默认值
                           tryCatch({
                             Model_IVW_BtoA_NoFilter_Test <- lm(Beta_A2[List_B] ~ -1 + Beta_B2[List_B])  # 使用所有显著SNP
                             P_Sim_IVW_BtoA_NoFilter_Test_k <- summary(Model_IVW_BtoA_NoFilter_Test)$coef[1,4]  # P值
                           }, error = function(e) {
                             P_Sim_IVW_BtoA_NoFilter_Test_k <- 1  # 出错时设为1
                           })
                           
                           # 测试集IVW分析（B→A方向，删除多效性SNP）
                           P_Sim_IVW_BtoA_Filter_Test_k <- 1  # 默认值
                           tryCatch({
                             Model_IVW_BtoA_Filter_Test <- lm(Beta_A2[List_B & !List_A] ~ -1 + Beta_B2[List_B & !List_A])  # 排除同时显著的SNP
                             P_Sim_IVW_BtoA_Filter_Test_k <- summary(Model_IVW_BtoA_Filter_Test)$coef[1,4]  # P值
                           }, error = function(e) {
                             P_Sim_IVW_BtoA_Filter_Test_k <- 1  # 出错时设为1
                           })
                           
                           # 返回单次模拟的结果
                           c(P_Sim_AtoB_Overall_k, P_Sim_BtoA_Overall_k, 
                             P_Sim_IVW_AtoB_NoFilter_Test_k, P_Sim_IVW_AtoB_Filter_Test_k,
                             P_Sim_IVW_BtoA_NoFilter_Test_k, P_Sim_IVW_BtoA_Filter_Test_k,
                             PRS_IV_cor_A_k, PRS_IV_Pperms_A_k, PRS_IV_cor_B_k, PRS_IV_Pperms_B_k)
                         }
  
  # 将并行结果转换为向量
  P_Sim_AtoB_Overall <- sim_results[,1]
  P_Sim_BtoA_Overall <- sim_results[,2]
  P_Sim_IVW_AtoB_NoFilter_Test <- sim_results[,3]
  P_Sim_IVW_AtoB_Filter_Test <- sim_results[,4]
  P_Sim_IVW_BtoA_NoFilter_Test <- sim_results[,5]
  P_Sim_IVW_BtoA_Filter_Test <- sim_results[,6]
  PRS_IV_cor_A <- sim_results[,7]
  PRS_IV_Pperms_A <- sim_results[,8]
  PRS_IV_cor_B <- sim_results[,9]
  PRS_IV_Pperms_B <- sim_results[,10]
  
  # 保存训练集数据
  train_data <- list(
    Data = Data,
    PhenoA = PhenoA,
    PhenoB = PhenoB,
    Beta_A = Beta_A,
    Beta_B = Beta_B,
    P_A = P_A,
    P_B = P_B,
    List_A = List_A,
    List_B = List_B,
    SNP_Effect_A = SNP_Effect_A,
    SNP_Effect_B = SNP_Effect_B,
    SNP_Pleiotropy_Effect_A = SNP_Pleiotropy_Effect_A,
    SNP_Pleiotropy_Effect_B = SNP_Pleiotropy_Effect_B,
    Keep_A = Keep_A,
    Keep_B = Keep_B,
    IVW_AtoB_NoFilter_Estimate_Train = IVW_AtoB_NoFilter_Estimate_Train,
    IVW_AtoB_Filter_Estimate_Train = IVW_AtoB_Filter_Estimate_Train,
    IVW_AtoB_NoFilter_P_Train = IVW_AtoB_NoFilter_P_Train,
    IVW_AtoB_Filter_P_Train = IVW_AtoB_Filter_P_Train,
    IVW_BtoA_NoFilter_Estimate_Train = IVW_BtoA_NoFilter_Estimate_Train,
    IVW_BtoA_Filter_Estimate_Train = IVW_BtoA_Filter_Estimate_Train,
    IVW_BtoA_NoFilter_P_Train = IVW_BtoA_NoFilter_P_Train,
    IVW_BtoA_Filter_P_Train = IVW_BtoA_Filter_P_Train,
    N_Ind_Train = N_Ind_Train,
    N_Ind_Test = N_Ind_Test,
    N_SNP = N_SNP,
    AtoB = AtoB,
    Pleio_Effect = pleio_effect,
    N_Eff_SNP_PhenoAB = N_Eff_SNP_PhenoAB
  )
  
  # 为训练集生成唯一文件名并保存到参数子文件夹
  train_filename <- paste0("train_data_N", N_Ind_Train, "_AtoB", AtoB, "_Pleio", pleio_effect, "_SNP", N_Eff_SNP_PhenoAB, ".RData")
  train_filepath <- file.path(folder_path, train_filename)
  save(train_data, file = train_filepath)
  
  # 保存测试集数据
  test_data <- list(
    P_Sim_AtoB_Overall = P_Sim_AtoB_Overall,
    P_Sim_BtoA_Overall = P_Sim_BtoA_Overall,
    P_Sim_IVW_AtoB_NoFilter_Test = P_Sim_IVW_AtoB_NoFilter_Test,
    P_Sim_IVW_AtoB_Filter_Test = P_Sim_IVW_AtoB_Filter_Test,
    P_Sim_IVW_BtoA_NoFilter_Test = P_Sim_IVW_BtoA_NoFilter_Test,
    P_Sim_IVW_BtoA_Filter_Test = P_Sim_IVW_BtoA_Filter_Test,
    PRS_IV_cor_A = PRS_IV_cor_A,
    PRS_IV_Pperms_A = PRS_IV_Pperms_A,
    PRS_IV_cor_B = PRS_IV_cor_B,
    PRS_IV_Pperms_B = PRS_IV_Pperms_B
  )
  
  # 返回结果
  return(list(
    # A→B方向结果
    PRS_AtoB_Power = mean(P_Sim_AtoB_Overall < 0.05),  # PRS A→B统计效力
    IVW_AtoB_NoFilter_Power_Test = mean(P_Sim_IVW_AtoB_NoFilter_Test < 0.05),  # IVW A→B不删SNP统计效力
    IVW_AtoB_Filter_Power_Test = mean(P_Sim_IVW_AtoB_Filter_Test < 0.05),  # IVW A→B删SNP统计效力
    IVW_AtoB_NoFilter_Estimate_Train = IVW_AtoB_NoFilter_Estimate_Train,  # IVW A→B不删SNP因果效应估计
    IVW_AtoB_Filter_Estimate_Train = IVW_AtoB_Filter_Estimate_Train,  # IVW A→B删SNP因果效应估计
    IVW_AtoB_NoFilter_P_Train = IVW_AtoB_NoFilter_P_Train,  # IVW A→B不删SNP P值
    IVW_AtoB_Filter_P_Train = IVW_AtoB_Filter_P_Train,  # IVW A→B删SNP P值
    
    # B→A方向结果
    PRS_BtoA_Power = mean(P_Sim_BtoA_Overall < 0.05),  # PRS B→A统计效力
    IVW_BtoA_NoFilter_Power_Test = mean(P_Sim_IVW_BtoA_NoFilter_Test < 0.05),  # IVW B→A不删SNP统计效力
    IVW_BtoA_Filter_Power_Test = mean(P_Sim_IVW_BtoA_Filter_Test < 0.05),  # IVW B→A删SNP统计效力
    IVW_BtoA_NoFilter_Estimate_Train = IVW_BtoA_NoFilter_Estimate_Train,  # IVW B→A不删SNP因果效应估计
    IVW_BtoA_Filter_Estimate_Train = IVW_BtoA_Filter_Estimate_Train,  # IVW B→A删SNP因果效应估计
    IVW_BtoA_NoFilter_P_Train = IVW_BtoA_NoFilter_P_Train,  # IVW B→A不删SNP P值
    IVW_BtoA_Filter_P_Train = IVW_BtoA_Filter_P_Train,  # IVW B→A删SNP P值
    
    # PRS工具变量验证结果（A→B方向）
    PRS_IV_cor_A_Mean = mean(PRS_IV_cor_A),  # PRS A与表型A平均相关性
    PRS_IV_Pperms_A_Mean = mean(PRS_IV_Pperms_A),  # PRS A工具变量permutation平均P值
    PRS_IV_Validation_Pass_Rate_A = mean(PRS_IV_Pperms_A < 0.05),  # PRS A工具变量验证通过率
    
    # PRS工具变量验证结果（B→A方向）
    PRS_IV_cor_B_Mean = mean(PRS_IV_cor_B),  # PRS B与表型B平均相关性
    PRS_IV_Pperms_B_Mean = mean(PRS_IV_Pperms_B),  # PRS B工具变量permutation平均P值
    PRS_IV_Validation_Pass_Rate_B = mean(PRS_IV_Pperms_B < 0.05),  # PRS B工具变量验证通过率
    
    # 训练集和测试集数据
    train_data = train_data,
    test_data = test_data
  ))
}

# ===================== Complete Batch Sensitivity Analysis =====================
# 参数设置
test_samples <- c(1000, 2000, 5000, 10000)  # 测试样本数量
AtoB_values <- c(0, 0.1, 0.2, 0.3)  # 因果效应大小
pleio_effects <- c(0.5, 1, 1.5, 2)  # 多效性效应大小
pleio_snp_counts <- c(10, 20, 30, 40, 50)  # 多效性SNP数量

# 初始化结果存储
results <- data.frame()  # 结果数据框

# 初始化数据存储列表
all_train_data <- list()
all_test_data <- list()
all_detailed_results <- list()

# 计算总组合数
total_combinations <- length(test_samples) * length(AtoB_values) * length(pleio_effects) * length(pleio_snp_counts)
current_combination <- 0  # 当前组合计数

cat("Starting complete sensitivity analysis with", total_combinations, "combinations...\n")

# 开始批处理分析
for (N_Ind_Test in test_samples) {  # 遍历测试样本数量
  for (AtoB in AtoB_values) {  # 遍历因果效应大小
    for (pleio_effect in pleio_effects) {  # 遍历多效性效应大小
      for (N_Eff_SNP_PhenoAB in pleio_snp_counts) {  # 遍历多效性SNP数量
        current_combination <- current_combination + 1  # 更新组合计数
        N_Ind_Train <- 20000  # 固定训练样本数量
        
        # 显示进度信息
        cat("Progress:", current_combination, "/", total_combinations, 
            " - N_Test=", N_Ind_Test, "AtoB=", AtoB, 
            "Pleio_Effect=", pleio_effect, "Pleio_SNP=", N_Eff_SNP_PhenoAB, "\n")
        
        # 创建参数对应的文件夹
        folder_info <- create_parameter_folder(N_Ind_Train, N_Ind_Test, AtoB, pleio_effect, N_Eff_SNP_PhenoAB)
        folder_path <- folder_info$main_folder
        test_folder_path <- folder_info$test_folder
        
        # 运行模拟
        sim_result <- run_simulation(N_Ind_Train, N_Ind_Test, AtoB, pleio_effect, 
                                     N_Eff_SNP_PhenoAB, N_Sim=1000, seed=52, 
                                     folder_path=folder_path, test_folder_path=test_folder_path)
        
        # 生成唯一的标识符
        data_id <- paste0("N_Train_", N_Ind_Train, "_",
                          "N_Test_", N_Ind_Test, "_",
                          "AtoB_", AtoB, "_",
                          "Pleio_Effect_", pleio_effect, "_",
                          "Pleio_SNP_", N_Eff_SNP_PhenoAB)
        
        # 存储训练集数据
        all_train_data[[data_id]] <- sim_result$train_data
        
        # 存储测试集数据
        all_test_data[[data_id]] <- sim_result$test_data
        
        # 存储所有带Mean的详细结果
        detailed_results <- list(
          # 参数设置
          N_Train = N_Ind_Train,
          N_Test = N_Ind_Test,
          AtoB = AtoB,
          Pleio_Effect = pleio_effect,
          N_Pleio_SNP = N_Eff_SNP_PhenoAB,
          
          # 原始测试集数据向量（1000次模拟的原始结果）
          P_Sim_AtoB_Overall_Raw = sim_result$test_data$P_Sim_AtoB_Overall,
          P_Sim_BtoA_Overall_Raw = sim_result$test_data$P_Sim_BtoA_Overall,
          P_Sim_IVW_AtoB_NoFilter_Test_Raw = sim_result$test_data$P_Sim_IVW_AtoB_NoFilter_Test,
          P_Sim_IVW_AtoB_Filter_Test_Raw = sim_result$test_data$P_Sim_IVW_AtoB_Filter_Test,
          P_Sim_IVW_BtoA_NoFilter_Test_Raw = sim_result$test_data$P_Sim_IVW_BtoA_NoFilter_Test,
          P_Sim_IVW_BtoA_Filter_Test_Raw = sim_result$test_data$P_Sim_IVW_BtoA_Filter_Test,
          PRS_IV_cor_A_Raw = sim_result$test_data$PRS_IV_cor_A,
          PRS_IV_Pperms_A_Raw = sim_result$test_data$PRS_IV_Pperms_A,
          PRS_IV_cor_B_Raw = sim_result$test_data$PRS_IV_cor_B,
          PRS_IV_Pperms_B_Raw = sim_result$test_data$PRS_IV_Pperms_B,
          
          # 汇总统计结果（Mean值）
          PRS_AtoB_Power = sim_result$PRS_AtoB_Power,
          PRS_BtoA_Power = sim_result$PRS_BtoA_Power,
          IVW_AtoB_NoFilter_Power_Test = sim_result$IVW_AtoB_NoFilter_Power_Test,
          IVW_AtoB_Filter_Power_Test = sim_result$IVW_AtoB_Filter_Power_Test,
          IVW_BtoA_NoFilter_Power_Test = sim_result$IVW_BtoA_NoFilter_Power_Test,
          IVW_BtoA_Filter_Power_Test = sim_result$IVW_BtoA_Filter_Power_Test,
          PRS_IV_cor_A_Mean = sim_result$PRS_IV_cor_A_Mean,
          PRS_IV_Pperms_A_Mean = sim_result$PRS_IV_Pperms_A_Mean,
          PRS_IV_Validation_Pass_Rate_A = sim_result$PRS_IV_Validation_Pass_Rate_A,
          PRS_IV_cor_B_Mean = sim_result$PRS_IV_cor_B_Mean,
          PRS_IV_Pperms_B_Mean = sim_result$PRS_IV_Pperms_B_Mean,
          PRS_IV_Validation_Pass_Rate_B = sim_result$PRS_IV_Validation_Pass_Rate_B,
          
          # 训练集分析结果
          IVW_AtoB_NoFilter_Estimate_Train = sim_result$IVW_AtoB_NoFilter_Estimate_Train,
          IVW_AtoB_Filter_Estimate_Train = sim_result$IVW_AtoB_Filter_Estimate_Train,
          IVW_AtoB_NoFilter_P_Train = sim_result$IVW_AtoB_NoFilter_P_Train,
          IVW_AtoB_Filter_P_Train = sim_result$IVW_AtoB_Filter_P_Train,
          IVW_BtoA_NoFilter_Estimate_Train = sim_result$IVW_BtoA_NoFilter_Estimate_Train,
          IVW_BtoA_Filter_Estimate_Train = sim_result$IVW_BtoA_Filter_Estimate_Train,
          IVW_BtoA_NoFilter_P_Train = sim_result$IVW_BtoA_NoFilter_P_Train,
          IVW_BtoA_Filter_P_Train = sim_result$IVW_BtoA_Filter_P_Train,
          
          # 额外统计信息
          P_Sim_AtoB_Overall_Mean = mean(sim_result$test_data$P_Sim_AtoB_Overall),
          P_Sim_AtoB_Overall_SD = sd(sim_result$test_data$P_Sim_AtoB_Overall),
          P_Sim_BtoA_Overall_Mean = mean(sim_result$test_data$P_Sim_BtoA_Overall),
          P_Sim_BtoA_Overall_SD = sd(sim_result$test_data$P_Sim_BtoA_Overall),
          P_Sim_IVW_AtoB_NoFilter_Test_Mean = mean(sim_result$test_data$P_Sim_IVW_AtoB_NoFilter_Test),
          P_Sim_IVW_AtoB_NoFilter_Test_SD = sd(sim_result$test_data$P_Sim_IVW_AtoB_NoFilter_Test),
          P_Sim_IVW_AtoB_Filter_Test_Mean = mean(sim_result$test_data$P_Sim_IVW_AtoB_Filter_Test),
          P_Sim_IVW_AtoB_Filter_Test_SD = sd(sim_result$test_data$P_Sim_IVW_AtoB_Filter_Test),
          P_Sim_IVW_BtoA_NoFilter_Test_Mean = mean(sim_result$test_data$P_Sim_IVW_BtoA_NoFilter_Test),
          P_Sim_IVW_BtoA_NoFilter_Test_SD = sd(sim_result$test_data$P_Sim_IVW_BtoA_NoFilter_Test),
          P_Sim_IVW_BtoA_Filter_Test_Mean = mean(sim_result$test_data$P_Sim_IVW_BtoA_Filter_Test),
          P_Sim_IVW_BtoA_Filter_Test_SD = sd(sim_result$test_data$P_Sim_IVW_BtoA_Filter_Test),
          
          # PRS工具变量额外统计信息
          PRS_IV_cor_A_Mean_Calc = mean(sim_result$test_data$PRS_IV_cor_A),
          PRS_IV_cor_A_SD = sd(sim_result$test_data$PRS_IV_cor_A),
          PRS_IV_Pperms_A_Mean_Calc = mean(sim_result$test_data$PRS_IV_Pperms_A),
          PRS_IV_Pperms_A_SD = sd(sim_result$test_data$PRS_IV_Pperms_A),
          PRS_IV_cor_B_Mean_Calc = mean(sim_result$test_data$PRS_IV_cor_B),
          PRS_IV_cor_B_SD = sd(sim_result$test_data$PRS_IV_cor_B),
          PRS_IV_Pperms_B_Mean_Calc = mean(sim_result$test_data$PRS_IV_Pperms_B),
          PRS_IV_Pperms_B_SD = sd(sim_result$test_data$PRS_IV_Pperms_B)
        )
        
        all_detailed_results[[data_id]] <- detailed_results
        
        # 保存当前参数组合的详细结果到参数子文件夹
        detailed_filename <- paste0("detailed_results_N", N_Ind_Train, "_", N_Ind_Test, "_AtoB", AtoB, "_Pleio", pleio_effect, "_SNP", N_Eff_SNP_PhenoAB, ".RData")
        detailed_filepath <- file.path(folder_path, detailed_filename)
        save(detailed_results, file = detailed_filepath)
        
        # 存储结果到数据框
        results <- rbind(results, data.frame(
          N_Train = N_Ind_Train,  # 训练样本数量
          N_Test = N_Ind_Test,  # 测试样本数量
          AtoB = AtoB,  # 因果效应大小
          Pleio_Effect = pleio_effect,  # 多效性效应大小
          N_Pleio_SNP = N_Eff_SNP_PhenoAB,  # 多效性SNP数量
          
          
          # A→B方向结果
          PRS_AtoB_Power = sim_result$PRS_AtoB_Power,  # PRS A→B统计效力
          IVW_AtoB_NoFilter_Power_Test = sim_result$IVW_AtoB_NoFilter_Power_Test,  # IVW A→B不删SNP统计效力
          IVW_AtoB_Filter_Power_Test = sim_result$IVW_AtoB_Filter_Power_Test,  # IVW A→B删SNP统计效力
          IVW_AtoB_NoFilter_Estimate_Train = sim_result$IVW_AtoB_NoFilter_Estimate_Train,  # IVW A→B不删SNP因果效应估计
          IVW_AtoB_Filter_Estimate_Train = sim_result$IVW_AtoB_Filter_Estimate_Train,  # IVW A→B删SNP因果效应估计
          IVW_AtoB_NoFilter_P_Train = sim_result$IVW_AtoB_NoFilter_P_Train,  # IVW A→B不删SNP P值
          IVW_AtoB_Filter_P_Train = sim_result$IVW_AtoB_Filter_P_Train,  # IVW A→B删SNP P值
          
          # B→A方向结果
          PRS_BtoA_Power = sim_result$PRS_BtoA_Power,  # PRS B→A统计效力
          IVW_BtoA_NoFilter_Power_Test = sim_result$IVW_BtoA_NoFilter_Power_Test,  # IVW B→A不删SNP统计效力
          IVW_BtoA_Filter_Power_Test = sim_result$IVW_BtoA_Filter_Power_Test,  # IVW B→A删SNP统计效力
          IVW_BtoA_NoFilter_Estimate_Train = sim_result$IVW_BtoA_NoFilter_Estimate_Train,  # IVW B→A不删SNP因果效应估计
          IVW_BtoA_Filter_Estimate_Train = sim_result$IVW_BtoA_Filter_Estimate_Train,  # IVW B→A删SNP因果效应估计
          IVW_BtoA_NoFilter_P_Train = sim_result$IVW_BtoA_NoFilter_P_Train,  # IVW B→A不删SNP P值
          IVW_BtoA_Filter_P_Train = sim_result$IVW_BtoA_Filter_P_Train,  # IVW B→A删SNP P值
          
          # PRS工具变量验证结果（A→B方向）
          PRS_IV_cor_A_Mean = sim_result$PRS_IV_cor_A_Mean,  # PRS A与表型A平均相关性
          PRS_IV_Pperms_A_Mean = sim_result$PRS_IV_Pperms_A_Mean,  # PRS A工具变量permutation平均P值
          PRS_IV_Validation_Pass_Rate_A = sim_result$PRS_IV_Validation_Pass_Rate_A,  # PRS A工具变量验证通过率
          
          # PRS工具变量验证结果（B→A方向）
          PRS_IV_cor_B_Mean = sim_result$PRS_IV_cor_B_Mean,  # PRS B与表型B平均相关性
          PRS_IV_Pperms_B_Mean = sim_result$PRS_IV_Pperms_B_Mean,  # PRS B工具变量permutation平均P值
          PRS_IV_Validation_Pass_Rate_B = sim_result$PRS_IV_Validation_Pass_Rate_B  # PRS B工具变量验证通过率
        ))
        
        # 显示当前结果
        cat("Done: PRS_AtoB_Power=", sim_result$PRS_AtoB_Power, 
            "PRS_BtoA_Power=", sim_result$PRS_BtoA_Power,
            "IVW_AtoB_NoFilter_Power_Test=", sim_result$IVW_AtoB_NoFilter_Power_Test,
            "IVW_AtoB_Filter_Power_Test=", sim_result$IVW_AtoB_Filter_Power_Test,
            "IVW_BtoA_NoFilter_Power_Test=", sim_result$IVW_BtoA_NoFilter_Power_Test,
            "IVW_BtoA_Filter_Power_Test=", sim_result$IVW_BtoA_Filter_Power_Test,
            "PRS_IV_cor_A=", round(sim_result$PRS_IV_cor_A_Mean, 3),
            "PRS_IV_cor_B=", round(sim_result$PRS_IV_cor_B_Mean, 3),
            "PRS_IV_Validation_Pass_Rate_A=", sim_result$PRS_IV_Validation_Pass_Rate_A,
            "PRS_IV_Validation_Pass_Rate_B=", sim_result$PRS_IV_Validation_Pass_Rate_B,
            "\n\n")
        
        # 保存当前参数组合的汇总结果到CSV文件
        summary_filename <- paste0("summary_results_N", N_Ind_Train, "_", N_Ind_Test, "_AtoB", AtoB, "_Pleio", pleio_effect, "_SNP", N_Eff_SNP_PhenoAB, ".csv")
        summary_filepath <- file.path(folder_path, summary_filename)
        current_result <- data.frame(
          N_Train = N_Ind_Train,  # 训练样本数量
          N_Test = N_Ind_Test,  # 测试样本数量
          AtoB = AtoB,  # 因果效应大小
          Pleio_Effect = pleio_effect,  # 多效性效应大小
          N_Pleio_SNP = N_Eff_SNP_PhenoAB,  # 多效性SNP数量
          
          # A→B方向结果
          PRS_AtoB_Power = sim_result$PRS_AtoB_Power,  # PRS A→B统计效力
          IVW_AtoB_NoFilter_Power_Test = sim_result$IVW_AtoB_NoFilter_Power_Test,  # IVW A→B不删SNP统计效力
          IVW_AtoB_Filter_Power_Test = sim_result$IVW_AtoB_Filter_Power_Test,  # IVW A→B删SNP统计效力
          IVW_AtoB_NoFilter_Estimate_Train = sim_result$IVW_AtoB_NoFilter_Estimate_Train,  # IVW A→B不删SNP因果效应估计
          IVW_AtoB_Filter_Estimate_Train = sim_result$IVW_AtoB_Filter_Estimate_Train,  # IVW A→B删SNP因果效应估计
          IVW_AtoB_NoFilter_P_Train = sim_result$IVW_AtoB_NoFilter_P_Train,  # IVW A→B不删SNP P值
          IVW_AtoB_Filter_P_Train = sim_result$IVW_AtoB_Filter_P_Train,  # IVW A→B删SNP P值
          
          # B→A方向结果
          PRS_BtoA_Power = sim_result$PRS_BtoA_Power,  # PRS B→A统计效力
          IVW_BtoA_NoFilter_Power_Test = sim_result$IVW_BtoA_NoFilter_Power_Test,  # IVW B→A不删SNP统计效力
          IVW_BtoA_Filter_Power_Test = sim_result$IVW_BtoA_Filter_Power_Test,  # IVW B→A删SNP统计效力
          IVW_BtoA_NoFilter_Estimate_Train = sim_result$IVW_BtoA_NoFilter_Estimate_Train,  # IVW B→A不删SNP因果效应估计
          IVW_BtoA_Filter_Estimate_Train = sim_result$IVW_BtoA_Filter_Estimate_Train,  # IVW B→A删SNP因果效应估计
          IVW_BtoA_NoFilter_P_Train = sim_result$IVW_BtoA_NoFilter_P_Train,  # IVW B→A不删SNP P值
          IVW_BtoA_Filter_P_Train = sim_result$IVW_BtoA_Filter_P_Train,  # IVW B→A删SNP P值
          
          # PRS工具变量验证结果（A→B方向）
          PRS_IV_cor_A_Mean = sim_result$PRS_IV_cor_A_Mean,  # PRS A与表型A平均相关性
          PRS_IV_Pperms_A_Mean = sim_result$PRS_IV_Pperms_A_Mean,  # PRS A工具变量permutation平均P值
          PRS_IV_Validation_Pass_Rate_A = sim_result$PRS_IV_Validation_Pass_Rate_A,  # PRS A工具变量验证通过率
          
          # PRS工具变量验证结果（B→A方向）
          PRS_IV_cor_B_Mean = sim_result$PRS_IV_cor_B_Mean,  # PRS B与表型B平均相关性
          PRS_IV_Pperms_B_Mean = sim_result$PRS_IV_Pperms_B_Mean,  # PRS B工具变量permutation平均P值
          PRS_IV_Validation_Pass_Rate_B = sim_result$PRS_IV_Validation_Pass_Rate_B  # PRS B工具变量验证通过率
        )
        write.csv(current_result, summary_filepath, row.names = FALSE)
        
        # 定期保存中间结果
        if (current_combination %% 10 == 0) {
          write.csv(results, "statistical_power_AposBneg1000_temp.csv", row.names = FALSE)
          cat("Intermediate results saved at combination", current_combination, "\n")
        }
      }
    }
  }
}

# 保存最终结果
write.csv(results, "statistical_power_AposBneg1000.csv", row.names = FALSE)
cat("\nAll complete sensitivity analysis finished! Results saved to statistical_power_AposBneg1000.csv\n")

# 保存所有数据到三个RData文件
cat("\nSaving all data to RData files...\n")

# 保存所有训练集数据
save(all_train_data, file = "all_train_data1000.RData")
cat("All training data saved to: all_train_data1000.RData\n")

# 保存所有测试集数据
save(all_test_data, file = "all_test_data1000.RData")
cat("All test data saved to: all_test_data1000.RData\n")

# 保存所有详细结果
save(all_detailed_results, file = "all_detailed_results1000.RData")
cat("All detailed results saved to: all_detailed_results1000.RData\n")

# 输出结果摘要
cat("\n=== Results Summary ===\n")
cat("Total simulations:", nrow(results), "\n")
cat("PRS A→B Power range:", range(results$PRS_AtoB_Power), "\n")
cat("PRS B→A Power range:", range(results$PRS_BtoA_Power), "\n")
cat("IVW A→B NoFilter Power (Test) range:", range(results$IVW_AtoB_NoFilter_Power_Test), "\n")
cat("IVW A→B Filter Power (Test) range:", range(results$IVW_AtoB_Filter_Power_Test), "\n")
cat("IVW B→A NoFilter Power (Test) range:", range(results$IVW_BtoA_NoFilter_Power_Test), "\n")
cat("IVW B→A Filter Power (Test) range:", range(results$IVW_BtoA_Filter_Power_Test), "\n")

# 按多效性效应分组显示结果
cat("\n=== Results by Pleiotropy Effect ===\n")
for (pleio in pleio_effects) {
  subset_results <- results[results$Pleio_Effect == pleio, ]
  cat("Pleio Effect =", pleio, ":\n")
  cat("  PRS A→B Power (mean):", mean(subset_results$PRS_AtoB_Power), "\n")
  cat("  PRS B→A Power (mean):", mean(subset_results$PRS_BtoA_Power), "\n")
  cat("  IVW A→B NoFilter Power (Test, mean):", mean(subset_results$IVW_AtoB_NoFilter_Power_Test), "\n")
  cat("  IVW A→B Filter Power (Test, mean):", mean(subset_results$IVW_AtoB_Filter_Power_Test), "\n")
  cat("  IVW B→A NoFilter Power (Test, mean):", mean(subset_results$IVW_BtoA_NoFilter_Power_Test), "\n")
  cat("  IVW B→A Filter Power (Test, mean):", mean(subset_results$IVW_BtoA_Filter_Power_Test), "\n")
}

# 按多效性SNP数量分组显示结果
cat("\n=== Results by Pleiotropy SNP Count ===\n")
for (pleio_count in pleio_snp_counts) {
  subset_results <- results[results$N_Pleio_SNP == pleio_count, ]
  cat("Pleio SNP Count =", pleio_count, ":\n")
  cat("  PRS A→B Power (mean):", mean(subset_results$PRS_AtoB_Power), "\n")
  cat("  PRS B→A Power (mean):", mean(subset_results$PRS_BtoA_Power), "\n")
  cat("  IVW A→B NoFilter Power (Test, mean):", mean(subset_results$IVW_AtoB_NoFilter_Power_Test), "\n")
  cat("  IVW A→B Filter Power (Test, mean):", mean(subset_results$IVW_AtoB_Filter_Power_Test), "\n")
  cat("  IVW B→A NoFilter Power (Test, mean):", mean(subset_results$IVW_BtoA_NoFilter_Power_Test), "\n")
  cat("  IVW B→A Filter Power (Test, mean):", mean(subset_results$IVW_BtoA_Filter_Power_Test), "\n")
}

# 按测试样本数量分组显示结果
cat("\n=== Results by Test Sample Size ===\n")
for (test_size in test_samples) {
  subset_results <- results[results$N_Test == test_size, ]
  cat("Test Sample Size =", test_size, ":\n")
  cat("  PRS A→B Power (mean):", mean(subset_results$PRS_AtoB_Power), "\n")
  cat("  PRS B→A Power (mean):", mean(subset_results$PRS_BtoA_Power), "\n")
  cat("  IVW A→B NoFilter Power (Test, mean):", mean(subset_results$IVW_AtoB_NoFilter_Power_Test), "\n")
  cat("  IVW A→B Filter Power (Test, mean):", mean(subset_results$IVW_AtoB_Filter_Power_Test), "\n")
  cat("  IVW B→A NoFilter Power (Test, mean):", mean(subset_results$IVW_BtoA_NoFilter_Power_Test), "\n")
  cat("  IVW B→A Filter Power (Test, mean):", mean(subset_results$IVW_BtoA_Filter_Power_Test), "\n")
}


# 显示保存文件总结
cat("\n=== File Saving Summary ===\n")
cat("Main results CSV file: statistical_power_AposBneg1000.csv\n")
cat("All training data saved in: all_train_data1000.RData\n")
cat("All test data saved in: all_test_data1000.RData\n")
cat("All detailed results saved in: all_detailed_results1000.RData\n")
cat("\nFolder structure created:\n")
cat("- Parameter folders: ", total_combinations, "folders\n")
cat("  Each folder contains:\n")
cat("    - train_data_N[params].RData (训练集数据)\n")
cat("    - detailed_results_N[params].RData (详细结果)\n")
cat("    - summary_results_N[params].csv (汇总结果)\n")
cat("    - test_sets/ (测试集子文件夹)\n")
cat("      - test_data_sim_[1-1000]_N[params].RData (1000个测试集)\n")
cat("\nTotal files created:", 4 + total_combinations * 3 + total_combinations * 1000, "files\n")
cat("- 1 main CSV file + 3 summary RData files +", total_combinations, "parameter folders\n")
cat("- Each parameter folder contains 3 files + 1000 test files\n")
cat("Folder naming convention: N_Train_[value]_N_Test_[value]_AtoB_[value]_Pleio_Effect_[value]_Pleio_SNP_[value]\n")
cat("File naming convention: [type]_N[params].RData/csv\n")
