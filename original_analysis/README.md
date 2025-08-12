# Original Analysis Workflow

This directory contains the original analysis scripts and workflows that were used to develop the Cox-MK R package.

## Contents

### Scripts
- `KnockoffScreen.r`: Original knockoff screening implementation
- `run_knockoff_gds.r`: Script for running knockoff analysis on GDS files

### Workflow Steps
- `step0 tte pheno/`: Time-to-event phenotype preparation
- `step1 Model fitting/`: Null model fitting with covariates  
- `step2 Knockoff generation/`: Generating knockoff variables
- `step3 SPA and variant association analysis/`: Association testing
- `step4 Calculate W statistics and select significant variants/`: Variable selection

### Results
- `plots/`: Generated plots and visualizations
- `simulation/`: Simulation studies and results

## Note

These files represent the original research workflow. The functionality has been refactored and packaged into the Cox-MK R package for easier use and distribution.

For current usage, please refer to the main Cox-MK package in the parent directory.
