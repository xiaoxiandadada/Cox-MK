# Cox-MK: Cox Regression with Model-X Knockoffs

A comprehensive R package for survival analysis using Cox regression with Model-X knockoffs for variable selection with false discovery rate (FDR) control.

## Features

- **Flexible Knockoff Generation**: User-controllable number of knockoff copies (M parameter)
- **SPACox Integration**: Optional support for large-scale genetic studies via [SPACox](https://github.com/WenjianBI/SPACox)
- **PLINK Data Support**: Direct loading of PLINK binary format data
- **FDR Control**: Rigorous variable selection with false discovery rate control
- **Survival Analysis**: Cox proportional hazards models for time-to-event data

## Installation

### From source (recommended)
```r
# Install dependencies
install.packages(c("Matrix", "survival", "irlba", "devtools"))

# Install CoxMK package
devtools::install_local("Cox-MK")
```

### From tar.gz
```r
install.packages("CoxMK_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick Start

```r
library(CoxMK)

# Load example data
data(example_genotypes)
data(example_positions) 
data(example_phenotype)
data(example_covariates)

# Create knockoffs (flexible M parameter)
knockoffs <- create_knockoffs(
  X = example_genotypes, 
  pos = example_positions, 
  M = 5  # or any positive integer
)

# Prepare phenotype data
pheno_data <- prepare_phenotype(example_phenotype, example_covariates)

# Fit null model (with optional SPACox)
null_model <- fit_null_model(
  time = pheno_data$time,
  status = pheno_data$status, 
  X = pheno_data[, c("age", "sex", "bmi")],
  use_spacox = TRUE  # Falls back to coxph if SPACox unavailable
)

# Run Cox knockoff screening
results <- cox_knockoff_screen(
  genotypes = example_genotypes,
  phenotype = pheno_data,
  null_model = null_model,
  knockoffs = knockoffs
)

# Calculate W statistics and select variables
w_stats <- calculate_w_statistics(results)
selected <- knockoff_filter(w_stats, fdr = 0.1)
```

## Project Structure

```
Cox-MK/
├── Cox-MK/              # R package source code
│   ├── R/               # Package functions
│   ├── data/            # Example datasets
│   ├── man/             # Documentation
│   ├── tests/           # Package tests
│   └── DESCRIPTION      # Package metadata
├── sample/              # Example PLINK data
├── original_analysis/   # Original analysis workflow
└── CoxMK_0.1.0.tar.gz  # Built package
```

## Key Functions

- `create_knockoffs()`: Generate Model-X knockoffs with flexible M parameter
- `fit_null_model()`: Fit Cox null model with SPACox support
- `cox_knockoff_screen()`: Perform association analysis
- `calculate_w_statistics()`: Compute W statistics for variable selection
- `knockoff_filter()`: Select variables with FDR control
- `load_plink_data()`: Load PLINK binary format data

## Citation

If you use this package in your research, please cite:

```bibtex
@software{CoxMK,
  title = {Cox-MK: Cox Regression with Model-X Knockoffs},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/xiaoxiandadada/Cox-MK}
}
```

## License

GPL-3

## References

- Candes, E., Fan, Y., Janson, L., & Lv, J. (2018). Panning for gold: Model-X knockoffs for high dimensional controlled variable selection. Journal of the Royal Statistical Society, 80(3), 551-577.
- Bi, W., et al. SPACox: Scalable and Accurate Cox Regression Analysis of Large-Scale Survival Data. https://github.com/WenjianBI/SPACox

