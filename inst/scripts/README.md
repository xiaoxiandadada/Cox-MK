# CoxMK User Script

This directory contains a user-friendly script for running Cox Model-X Knockoff analysis without needing to write R code.

## Script Available

### `run_coxmk.R` - Complete Analysis Script

A comprehensive script with command-line argument parsing using the `optparse` package.

**Installation Requirements:**
```r
install.packages("optparse")
```

**Usage:**
```bash
Rscript run_coxmk.R --plink_prefix data/sample --phenotype pheno.txt --covariates covar.txt [options]
```

**Required Arguments:**
- `--plink_prefix`: Path prefix for PLINK files (without .bed/.bim/.fam extension)
- `--phenotype`: Path to phenotype file with columns: IID, time, status
- `--covariates`: Path to covariates file with columns: IID, age, sex, etc.

**Optional Arguments:**
- `--M`: Number of knockoff copies (default: 5)
- `--fdr`: False discovery rate level (default: 0.05)
- `--method`: W statistics method: median, difference, ratio (default: median)
- `--output_dir`: Output directory for results (default: current directory)
- `--gds_file`: Pre-generated GDS file with knockoffs (optional)
- `--null_model`: Pre-fitted null model RDS file (optional)
- `--time_col`: Column name for survival time (default: time)
- `--status_col`: Column name for event status (default: status)
- `--help`: Show help message

**Examples:**

```bash
# 显示帮助信息
Rscript run_coxmk.R --help

# 基本分析（使用示例数据）
Rscript run_coxmk.R \
  --plink_prefix inst/extdata/sample \
  --phenotype inst/extdata/tte_phenotype.txt \
  --covariates inst/extdata/covariates.txt

# 自定义参数分析
Rscript run_coxmk.R \
  --plink_prefix data/ukb_sample \
  --phenotype phenotype.txt \
  --covariates covariates.txt \
  --M 10 \
  --fdr 0.1 \
  --method difference \
  --output_dir results/

# 快速分析（最少参数）
Rscript run_coxmk.R \
  --plink_prefix data/sample \
  --phenotype pheno.txt \
  --covariates covar.txt

# 使用GDS文件重用knockoffs
Rscript run_coxmk.R \
  --plink_prefix data/sample \
  --phenotype pheno.txt \
  --covariates covar.txt \
  --gds_file knockoffs.gds \
  --M 5 \
  --fdr 0.05
```

## Input File Formats

### Phenotype File
Tab-separated file with columns:
- `IID`: Individual ID (must match PLINK .fam file)
- `time`: Survival time (numeric, > 0)
- `status`: Event status (0 = censored, 1 = event)

Example:
```
IID	time	status
sample_1	10.5	1
sample_2	8.2	0
sample_3	15.7	1
```

### Covariates File
Tab-separated file with columns:
- `IID`: Individual ID (must match phenotype file)
- Additional columns: Numeric covariates (age, sex, PCs, etc.)

Example:
```
IID	age	sex	PC1	PC2
sample_1	45	1	0.12	-0.34
sample_2	52	0	-0.08	0.21
sample_3	38	1	0.15	0.09
```

### PLINK Files
Standard PLINK binary format:
- `.bed`: Binary genotype data
- `.bim`: Variant information
- `.fam`: Sample information

## Output Files

Both scripts generate the following output files:

1. **`coxmk_results.rds`** (or `results.rds` for quick script)
   - Complete R object with all analysis results
   - Can be loaded with `readRDS()` for further analysis

2. **`selected_snps.txt`** (if any SNPs are selected)
   - Tab-separated file with selected variants
   - Columns: index, snp_id, chromosome, position, W_statistic

3. **`analysis_summary.txt`**
   - Human-readable summary of the analysis
   - Input parameters, data summary, and results

4. **`knockoffs.gds`** (if GDS file option is used)
   - Binary file containing generated knockoff matrices
   - Can be reused for subsequent analyses

## Workflow Examples

### Basic Workflow
```bash
# 1. Prepare your data files
# 2. Run analysis
Rscript run_coxmk.R --plink_prefix data/genotypes --phenotype survival_data.txt --covariates demographics.txt --output_dir results/

# 3. Check results
cat results/analysis_summary.txt
```

### Two-Step Workflow (with GDS file reuse)
```bash
# Step 1: Generate knockoffs and save for reuse
Rscript run_coxmk.R \
  --plink_prefix data/genotypes \
  --phenotype survival_data.txt \
  --covariates demographics.txt \
  --gds_file knockoffs.gds \
  --output_dir analysis1/

# Step 2: Reuse knockoffs for different phenotype/parameters
Rscript run_coxmk.R \
  --plink_prefix data/genotypes \
  --phenotype different_phenotype.txt \
  --covariates demographics.txt \
  --gds_file knockoffs.gds \
  --fdr 0.1 \
  --output_dir analysis2/
```

## Troubleshooting

### Common Issues

1. **"optparse package not found"**
   ```r
   install.packages("optparse")
   ```

2. **"PLINK files not found"**
   - Check that all three files (.bed, .bim, .fam) exist
   - Use absolute paths or ensure files are in the correct location

3. **"Sample ID mismatch"**
   - Ensure IID columns match between phenotype, covariates, and PLINK .fam files
   - Check for extra spaces or different ID formats

4. **Memory issues with large datasets**
   - Use smaller M (number of knockoffs)
   - Consider running on a machine with more RAM
   - Use the GDS workflow to save/reuse knockoffs

### Getting Help

For each script, you can get help by running:
```bash
# Full script
Rscript run_coxmk.R --help

# Quick script
Rscript quick_coxmk.R --help
```

## Performance Tips

1. **Start with small parameters**: Use M=5 and a subset of your data for initial testing
2. **Use GDS files**: Save knockoffs once and reuse for multiple analyses
3. **Monitor memory usage**: Large datasets may require more RAM
4. **Parallel processing**: The underlying functions use optimized linear algebra libraries
