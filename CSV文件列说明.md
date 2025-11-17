# pleiotropy_variance_explained_all_settings.csv 列说明

## 文件概览

- **总行数**：320行（4种设置 × 4个Pleio_Effect × 5个N_Pleio_SNP × 4个AtoB）
- **总列数**：17列

---

## 列说明

### 1. Setting（设置类型）
- **含义**：模拟设置类型
- **可能值**：
  - `BothPos`：多效性SNP对表型A和B均为正效应
  - `BothNeg`：多效性SNP对表型A和B均为负效应
  - `AposBneg`：多效性SNP对表型A为正效应，对表型B为负效应
  - `AnegBpos`：多效性SNP对表型A为负效应，对表型B为正效应

### 2. Pleio_Effect（多效性效应大小）
- **含义**：多效性SNP的效应大小参数
- **可能值**：0.5, 1, 1.5, 2
- **说明**：用于生成多效性SNP效应大小的均值（标准差等于均值）

### 3. N_Pleio_SNP（多效性SNP数量）
- **含义**：同时影响表型A和B的多效性SNP数量
- **可能值**：10, 20, 30, 40, 50

### 4. AtoB（因果效应系数）
- **含义**：表型A对表型B的因果效应系数
- **可能值**：0, 0.1, 0.2, 0.3
- **说明**：当AtoB = 0时，表示无因果效应；AtoB > 0表示存在正向因果效应

---

## 表型A的方差分解结果

### 5. PhenoA_Specific_Pct_Genetic（表型A - 特异性SNP占遗传方差百分比）
- **含义**：仅影响表型A的特异性SNP在表型A遗传方差中的占比
- **单位**：百分比（%）
- **范围**：通常为90-99%
- **计算公式**：`var(特异性SNP对表型A的贡献) / var(表型A的总遗传成分) × 100%`

### 6. PhenoA_Specific_Pct_Total（表型A - 特异性SNP占总表型方差百分比）
- **含义**：仅影响表型A的特异性SNP在表型A总方差（含环境噪声）中的占比
- **单位**：百分比（%）
- **范围**：通常为18-20%
- **计算公式**：`var(特异性SNP对表型A的贡献) / var(表型A的总方差) × 100%`

### 7. PhenoA_Pleio_Pct_Genetic（表型A - 多效性SNP占遗传方差百分比）
- **含义**：多效性SNP在表型A遗传方差中的占比
- **单位**：百分比（%）
- **范围**：通常为1-33%（取决于Pleio_Effect和N_Pleio_SNP）
- **计算公式**：`var(多效性SNP对表型A的贡献) / var(表型A的总遗传成分) × 100%`

### 8. PhenoA_Pleio_Pct_Total（表型A - 多效性SNP占总表型方差百分比）
- **含义**：多效性SNP在表型A总方差（含环境噪声）中的占比
- **单位**：百分比（%）
- **范围**：通常为0.1-0.6%
- **计算公式**：`var(多效性SNP对表型A的贡献) / var(表型A的总方差) × 100%`

### 9. PhenoA_Env_Pct_Total（表型A - 环境噪声占总表型方差百分比）
- **含义**：环境噪声在表型A总方差中的占比
- **单位**：百分比（%）
- **范围**：通常为79-80%
- **计算公式**：`var(环境噪声) / var(表型A的总方差) × 100%`
- **说明**：环境噪声方差 = 总表型方差 - 总遗传方差

---

## 表型B的方差分解结果

### 10. PhenoB_Causal_Pct_Genetic（表型B - A→B因果效应占遗传方差百分比）
- **含义**：A→B因果效应在表型B遗传方差中的占比
- **单位**：百分比（%）
- **范围**：
  - AtoB = 0时：0%
  - AtoB = 0.1时：约1%
  - AtoB = 0.2时：约3-4%
  - AtoB = 0.3时：约7-8%
- **计算公式**：`var(表型A的遗传成分 × AtoB) / var(表型B的总遗传成分) × 100%`

### 11. PhenoB_Causal_Pct_Total（表型B - A→B因果效应占总表型方差百分比）
- **含义**：A→B因果效应在表型B总方差（含环境噪声）中的占比
- **单位**：百分比（%）
- **范围**：通常为0-1.7%（取决于AtoB值）
- **计算公式**：`var(表型A的遗传成分 × AtoB) / var(表型B的总方差) × 100%`

### 12. PhenoB_Pleio_Pct_Genetic（表型B - 多效性SNP占遗传方差百分比）
- **含义**：多效性SNP在表型B遗传方差中的占比
- **单位**：百分比（%）
- **范围**：通常为1-37%（取决于Pleio_Effect和N_Pleio_SNP）
- **计算公式**：`var(多效性SNP对表型B的贡献) / var(表型B的总遗传成分) × 100%`

### 13. PhenoB_Pleio_Pct_Total（表型B - 多效性SNP占总表型方差百分比）
- **含义**：多效性SNP在表型B总方差（含环境噪声）中的占比
- **单位**：百分比（%）
- **范围**：通常为0.1-0.6%
- **计算公式**：`var(多效性SNP对表型B的贡献) / var(表型B的总方差) × 100%`

### 14. PhenoB_Specific_Pct_Genetic（表型B - 特异性SNP占遗传方差百分比）
- **含义**：仅影响表型B的特异性SNP在表型B遗传方差中的占比
- **单位**：百分比（%）
- **范围**：通常为77-99%（取决于其他成分的占比）
- **计算公式**：`var(特异性SNP对表型B的贡献) / var(表型B的总遗传成分) × 100%`

### 15. PhenoB_Specific_Pct_Total（表型B - 特异性SNP占总表型方差百分比）
- **含义**：仅影响表型B的特异性SNP在表型B总方差（含环境噪声）中的占比
- **单位**：百分比（%）
- **范围**：通常为18-20%
- **计算公式**：`var(特异性SNP对表型B的贡献) / var(表型B的总方差) × 100%`

### 16. PhenoB_Env_Pct_Total（表型B - 环境噪声占总表型方差百分比）
- **含义**：环境噪声在表型B总方差中的占比
- **单位**：百分比（%）
- **范围**：通常为79-80%
- **计算公式**：`var(环境噪声) / var(表型B的总方差) × 100%`
- **说明**：环境噪声方差 = 总表型方差 - 总遗传方差

---

## 数据完整性检查

### 百分比总和验证

**表型A（占遗传方差百分比）：**
- PhenoA_Specific_Pct_Genetic + PhenoA_Pleio_Pct_Genetic ≈ 100%
- **计算方法**：由于表型A的各成分（特异性SNP和多效性SNP）使用不同的SNP集合，它们之间相互独立，因此可以直接用方差比计算。

**表型A（占总表型方差百分比）：**
- PhenoA_Specific_Pct_Total + PhenoA_Pleio_Pct_Total + PhenoA_Env_Pct_Total ≈ 100%

**表型B（占遗传方差百分比）：**
- PhenoB_Causal_Pct_Genetic + PhenoB_Pleio_Pct_Genetic + PhenoB_Specific_Pct_Genetic = 100%
- **计算方法**：使用协方差方法计算各成分对总遗传方差的贡献：
  - `PhenoB_Causal_Pct_Genetic = cov(PhenoB_genetic, PhenoB_causal) / var(PhenoB_genetic) × 100%`
  - `PhenoB_Pleio_Pct_Genetic = cov(PhenoB_genetic, PhenoB_pleio) / var(PhenoB_genetic) × 100%`
  - `PhenoB_Specific_Pct_Genetic = cov(PhenoB_genetic, PhenoB_specific) / var(PhenoB_genetic) × 100%`
- **为什么需要协方差方法**：由于`PhenoB_causal = PhenoA_genetic × AtoB`，而`PhenoA_genetic`包含多效性SNP对A的效应，因此`PhenoB_causal`和`PhenoB_pleio`之间存在协方差。如果直接用各成分的方差除以总方差，总和可能不等于100%（可能小于或大于100%，取决于协方差的符号）。使用协方差方法可以确保总和严格等于100%。

**表型B（占总表型方差百分比）：**
- PhenoB_Causal_Pct_Total + PhenoB_Pleio_Pct_Total + PhenoB_Specific_Pct_Total + PhenoB_Env_Pct_Total ≈ 100%
- **计算方法**：同样使用协方差方法，但分母是总表型方差。

**验证结果**：
- 表型B占遗传方差百分比总和：均值100.00%，范围99.99%-100.01%
- 表型B占总表型方差百分比总和：均值99.97%，范围99.22%-100.71%

---

## 使用示例

### 查询特定参数组合的结果

例如，查询 BothPos 设置，Pleio_Effect = 1, N_Pleio_SNP = 30, AtoB = 0 的结果：

```r
df <- read.csv('pleiotropy_variance_explained_all_settings.csv')
result <- df[df$Setting == "BothPos" & 
            df$Pleio_Effect == 1 & 
            df$N_Pleio_SNP == 30 & 
            df$AtoB == 0, ]
```

结果将显示：
- 表型A：特异性SNP约92.4%，多效性SNP约7.6%（占遗传方差）
- 表型B：A→B因果效应0%，多效性SNP约7.9%，特异性SNP约92.0%（占遗传方差）

### 比较不同设置

```r
# 比较四种设置在相同参数下的结果
comparison <- df[df$Pleio_Effect == 1 & 
                 df$N_Pleio_SNP == 30 & 
                 df$AtoB == 0, ]
```

---

## 关键指标说明

### 占遗传方差百分比 vs 占总表型方差百分比

1. **占遗传方差百分比**：
   - 反映各遗传成分在遗传变异中的相对重要性
   - 总和应为100%（因为所有遗传成分之和等于总遗传成分）
   - 不受环境噪声影响

2. **占总表型方差百分比**：
   - 反映各成分在总表型变异中的贡献
   - 总和应为100%（遗传成分 + 环境噪声 = 总表型方差）
   - 受环境噪声影响（环境噪声通常占约80%）

### 为什么两种百分比都需要？

- **占遗传方差百分比**：用于理解遗传成分之间的相对重要性
- **占总表型方差百分比**：用于理解各成分在总表型变异中的实际贡献（更贴近真实世界的情况）

---

**文档生成日期**：2025年11月16日  
**数据文件**：pleiotropy_variance_explained_all_settings.csv

