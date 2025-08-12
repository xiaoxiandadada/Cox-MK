#!/bin/bash
#SBATCH --job-name=ss
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/step3_log/ss_%A_%a.out

############################ 用户可修改区域 ###################################
FILELIST="/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/gds_filelist.txt"
GDS_DIR=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/generate
RSCRIPT=/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/step3.R
OUTDIR=/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/step3_res_h
LOGDIR=/lustre/home/acct-mashiyang1991/dhhhh020606/tzh2/step3_log_h
STEP=500                     
###############################################################################

if [[ "$ROLE" == "worker" ]]; then
    module load miniconda3
    eval "$(conda shell.bash hook)"
    conda activate R  

    GDS_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST" | tr -d '\r')
    
    Rscript "$RSCRIPT" "$GDS_FILE" "$OUTDIR"
    
    exit 0
fi


total=$(wc -l < "$FILELIST")
for ((from=1; from<=total; from+=STEP)); do
    to=$(( from + STEP - 1 )); (( to > total )) && to=$total

    jobid=$(sbatch --parsable \
        --array=${from}-${to}%${STEP} \
        --export=ROLE=worker,FILELIST=$FILELIST,RSCRIPT=$RSCRIPT,OUTDIR=$OUTDIR \
        --output="$LOGDIR/slurm_%A_%a.out" \
        --error="$LOGDIR/slurm_%A_%a.err" \
        "$0")

    echo "Submitted batch $jobid for files $from-$to ($(date))"
    while squeue -h -j "$jobid" 2>/dev/null | grep -q . ; do sleep 60; done
done
echo "All batches finished ($(date))."