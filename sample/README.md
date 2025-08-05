# plink格式示例数据

## 文件说明
- **sample.bed**: plink二进制基因型文件 (25KB)
- **sample.bim**: SNP信息文件 (1000个变异位点)  
- **sample.fam**: 样本信息文件 (100个样本)
- **tte_phenotype.txt**: Time-to-event表型数据
- **covariates.txt**: 协变量数据 (年龄、性别、BMI、吸烟)

## 数据特征
- 样本数: 100
- SNP数: 1000 (分布在1-22号染色体)
- 事件发生率: ~31%
- MAF范围: 0.05-0.5

## 使用方法
```bash
# plink质控和分析
plink --bfile sample --freq --out sample_freq
plink --bfile sample --maf 0.05 --make-bed --out sample_qc
```

```r
# R中读取
library(data.table)
tte <- fread("tte_phenotype.txt")
cov <- fread("covariates.txt")
```
