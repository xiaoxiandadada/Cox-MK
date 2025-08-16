<div align="center">
  <img src="archive/coxmk_capsule_strong.svg" alt="CoxMK" width="80%"/>
</div>

[![R](https://img.shields.io/badge/R-%3E%3D3.5.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

CoxMK implements Cox regression with Multiple knockoffs for survival analysis, providing finite-sample FDR control in high-dimensional genetic studies. 

**Features:**
- 🧬 PLINK format support for large-scale GWAS
- 📊 Finite-sample FDR control via Mutiple knockoffs  
- ⚡ Efficient null model fitting using **SPACox** for large datasets
- 🖥️ Command line interface for easy use

## Installation

```r
# Install from GitHub
devtools::install_github("xiaoxiandadada/Cox-MK")

# Install dependencies
install.packages(c("Matrix", "survival", "irlba", "optparse"))

# Install SPACox for efficient null model fitting (recommended)
devtools::install_github("WenjianBI/SPACox")
```

## Quick Start

### R Interface

```r
library(CoxMK)

# Load example data
extdata_path <- system.file("extdata", package = "CoxMK")

# Run analysis
result <- cox_knockoff_analysis(
  plink_prefix = file.path(extdata_path, "sample"),
  phenotype_file = file.path(extdata_path, "tte_phenotype.txt"),
  covariate_file = file.path(extdata_path, "covariates.txt"),
  M = 5,           # Number of knockoffs
  fdr = 0.05       # False discovery rate
)

print(result$summary)
```

### Command Line Interface

```bash
# Show help
Rscript inst/scripts/run_coxmk.R --help

# Run analysis
Rscript inst/scripts/run_coxmk.R \
  --plink_prefix data/sample \
  --phenotype phenotype.txt \
  --covariates covariates.txt \
  --M 5 --fdr 0.05 --output_dir results/
```

## Input Data

**Required files:**
- PLINK binary files (`.bed`, `.bim`, `.fam`)
- Phenotype file: `IID`, `time`, `status` columns
- Covariates file: `IID` + covariate columns

**Example phenotype file:**
```
IID      time    status
sample1  10.5    1
sample2  8.2     0
```

## Output

```r
result$selected_vars   # Selected SNP indices
result$W_stats        # W statistics for all SNPs
result$summary        # Analysis summary
```

## Citation

```
@software{CoxMK2025,
  title   = {CoxMK: Cox Regression with Mutiple Knockoffs for Survival Analysis},
  author  = {Chen, Yang},
  year    = {2025},
  version = {0.1.0},
  url     = {https://github.com/xiaoxiandadada/Cox-MK},
  note    = {R package}
}
```

## Contact

For questions, please contact yangchen5@stu.scu.edu.cn.

## License

This software is licensed under GPLv3.
