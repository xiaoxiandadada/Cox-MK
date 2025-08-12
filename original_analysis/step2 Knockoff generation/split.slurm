#!/bin/bash
#SBATCH --job-name=cut_plink
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/log/cut_%A_%a.out


# ------------- 路径 -------------
COORD_FILE="/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/coords.txt"
IN_PREFIX="/lustre/home/acct-mashiyang1991/dhhhh020606/data/imputed_data_QC/ukb_imputed_QC_wb_white_unrelated_maf001"
OUT_DIR="/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/split"
LOG_DIR="/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/log"


# ------------- 读取坐标 -------------
read CHR START STOP <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" "$COORD_FILE")
REGION_TAG="${CHR}_${START}_${STOP}"
OUT_PREFIX="${OUT_DIR}/ukb_split_${REGION_TAG}"


plink2 \
  --bfile "$IN_PREFIX" \
  --chr "$CHR" \
  --from-bp "$START" \
  --to-bp   "$STOP" \
  --make-bed \
  --out "$OUT_PREFIX" \
  --threads $SLURM_CPUS_PER_TASK

echo "[$(date)] Finished cutting region ${REGION_TAG}"