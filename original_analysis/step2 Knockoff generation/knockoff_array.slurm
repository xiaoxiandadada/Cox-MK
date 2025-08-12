#!/bin/bash
#SBATCH --job-name=knockoff
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/genelog/ko_%A_%a.out

############################ 用户可修改区域 ###################################
COORDS=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/coords.txt
RSCRIPT=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/run_knockoff_gds.R
OUTDIR=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/generate
LOGDIR=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/genelog
STEP=500                       # 每批最大并发

###############################################################################

# ───────────────────────── Worker 功能 ───────────────────────────────────────
if [[ "$ROLE" == "worker" ]]; then
    module load miniconda3
    eval "$(conda shell.bash hook)"
    conda activate R

    # 提取并清理一行（去掉潜在的 ^M）
    line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BLOCKFILE" | tr -d '\r')
    CHR=$(awk '{print $1}' <<<"$line")
    START=$(awk '{print $2}' <<<"$line")
    END=$(awk '{print $3}' <<<"$line")

    Rscript "$RSCRIPT" "$CHR" "$START" "$END"
    exit 0
fi




total=$(wc -l < "$COORDS")                              # 区块总数 2495
for ((from=1; from<=total; from+=STEP)); do
    to=$(( from + STEP - 1 )); (( to > total )) && to=$total

    jobid=$(sbatch --parsable \
        --array=${from}-${to}%${STEP} \
        --export=ROLE=worker,BLOCKFILE=$COORDS,RSCRIPT=$RSCRIPT \
        "$0")

    echo "Submitted batch $jobid for lines $from-$to ($(date))"

    # 等此批所有子任务完成，再提交下一批
    while squeue -h -j "$jobid" 2>/dev/null | grep -q . ; do sleep 60; done
done
echo "All batches finished ($(date))."
