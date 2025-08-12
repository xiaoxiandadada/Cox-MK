#!/bin/bash
#SBATCH --job-name=fl
#SBATCH --partition=192c6t
#SBATCH --mem=2T
#SBATCH --ntasks-per-node=10
#SBATCH --output=slurm-%j.out

plink2 \
  --bed /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb.bed \
  --bim /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb.bim \
  --fam /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb.fam \
  --keep /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/white_unrelated_ids_plink.txt \
   --maf 0.01 \
   --make-bed \
   --out /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb_white_unrelated_maf001
    --freq \
  --out /lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb_white_unrelated_freq