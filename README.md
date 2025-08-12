# Cox-MK

Cox Regression with Model-X Knockoffs for Survival Analysis

## Overview

`CoxMK` implements Cox regression analysis with Model-X knockoffs for variable selection in survival analysis. The package provides functions for:

- Generating knockoff variables for genetic data (user-configurable number of knockoffs)
- Fitting Cox proportional hazards models (with optional SPACox support)
- Performing association analysis with FDR control
- Variable selection using the knockoff filter

## Key Features

### SPACox Integration
For large-scale genetic studies, this package supports [SPACox](https://github.com/WenjianBI/SPACox) for efficient null model fitting. SPACox is particularly recommended for genome-wide association studies with survival outcomes.

### Flexible Knockoff Generation
Users can specify the number of knockoff copies (M) based on their computational resources and power requirements:

## Installation

```r
# Install development version from GitHub
devtools::install_github("xiaoxiandadada/Cox-MK")

# Optional: Install SPACox for large-scale studies
devtools::install_github("WenjianBI/SPACox")
```

## Quick Example

```r
library(CoxMK)

# Load example data
data(example_genotypes)
data(example_positions)
data(example_phenotype)
data(example_covariates)

# Create knockoff variables (user can choose M)
knockoffs_5 <- create_knockoffs(
  X = example_genotypes,
  pos = example_positions,
  M = 5  # Default: 5 knockoffs, users can specify 3, 10, etc.
)

# Alternative: fewer knockoffs for faster computation
knockoffs_3 <- create_knockoffs(
  X = example_genotypes,
  pos = example_positions,
  M = 3  # Faster but potentially less powerful
)

# Prepare phenotype data
pheno_data <- merge(example_phenotype, example_covariates, by = c("FID", "IID"))
covariates <- pheno_data[, c("age", "sex", "bmi")]

# Fit null model with standard coxph
null_model <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates
)

# For large-scale studies, use SPACox (if installed)
null_model_spa <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status,
  covariates = covariates,
  use_spacox = TRUE  # Requires SPACox package
)

# Perform knockoff screening (for demonstration, using simulated results)
set.seed(123)
p <- ncol(example_genotypes)
original_pvals <- runif(p, 0.001, 1)
original_coefs <- rnorm(p, 0, 0.5)
knockoff_pvals <- lapply(1:5, function(k) runif(p, 0.001, 1))  # 5 knockoffs
knockoff_coefs <- lapply(1:5, function(k) rnorm(p, 0, 0.5))

# Calculate W statistics and apply filter
w_stats <- calculate_w_statistics(original_pvals, knockoff_pvals, 
                                 original_coefs, knockoff_coefs)
selected <- knockoff_filter(w_stats, fdr = 0.1)

print(paste("Selected", length(selected$selected), "variables"))
```

## Features

### Data Loading
- `load_plink_data()`: Load PLINK binary format data
- `prepare_phenotype()`: Prepare time-to-event phenotype data

### Knockoff Generation
- `create_knockoffs()`: Generate Model-X knockoff variables with:
  - Leveraging scores for importance sampling
  - LD-based clustering
  - Regression-based knockoff construction

### Cox Regression Analysis
- `fit_null_model()`: Fit null Cox model with covariates
- `cox_knockoff_screen()`: Perform association analysis
- `calculate_w_statistics()`: Calculate W statistics for selection
- `knockoff_filter()`: Apply knockoff filter with FDR control

## Input Data Format

### Genotype Data
- PLINK binary format (.bed/.bim/.fam) or sparse matrices
- Samples in rows, SNPs in columns
- 0/1/2 encoding for genotype counts

### Phenotype Data
Text file with columns:
- `FID`: Family ID
- `IID`: Individual ID  
- `time`: Survival time
- `status`: Event indicator (0/1)

### Covariates
Text file with columns:
- `FID`, `IID`: Sample identifiers
- Additional covariate columns (age, sex, etc.)

## Method Overview

The knockoff method provides finite-sample FDR control by:

1. **Generating knockoffs**: Create synthetic variables that mimic the correlation structure of original variables
2. **Variable importance**: Test both original and knockoff variables
3. **W statistics**: Calculate signed maximum or difference statistics
4. **Selection**: Apply data-adaptive threshold to control FDR

## Citation

If you use this package, please cite:

```
CoxMK: Cox Regression with Model-X Knockoffs for Survival Analysis
Yang Chen (yangchen5&#64;stu.scu.edu.cn)
R package version 0.1.0 (2025)
https://github.com/xiaoxiandadada/Cox-MK
```

Or in BibTeX format:

```bibtex
@software{CoxMK2025,
  title = {CoxMK: Cox Regression with Model-X Knockoffs for Survival Analysis},
  author = {Yang Chen},
  year = {2025},
  url = {https://github.com/xiaoxiandadada/Cox-MK},
  note = {R package version 0.1.0}
}
```

## Contact

For questions or issues, please contact Yang Chen at yangchen5&#64;stu.scu.edu.cn

## License

GPL-3

