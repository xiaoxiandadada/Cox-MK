#!/bin/bash
#SBATCH --job-name=null
#SBATCH --partition=debug
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1


module load miniconda3
source activate R

data_dir=/lustre/home/acct-mashiyang1991/dhhhh020606/hl/Cox-knockoff/code/step1/data

Rscript fit_null.R \
  "${data_dir}/pheno_white_unrelated.csv" \
  "${data_dir}/covar_white_unrelated.txt" \
  "${data_dir}/white_unrelated_ids.txt" \
  "${data_dir}/ukb_imputed_QC_wb_white_unrelated_maf001.fam" \
  //lustre/home/acct-mashiyang1991/dhhhh020606/hl/Cox-knockoff/code/step1/data