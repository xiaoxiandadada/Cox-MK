#!/bin/bash
#SBATCH --job-name=master_cut
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/log/master_cut_%j.out
#SBATCH --error=/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/log/master_cut_%j.err


# 参数
step=500
total=2495
coords="/lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/coords.txt"

for ((start=1; start<=total; start+=step)); do
    end=$((start+step-1))
    [ $end -gt $total ] && end=$total
    echo "[$(date)] Submitting array ${start}-${end}"
    jobid=$(sbatch --array=${start}-${end} \
                   --parsable \
                   /lustre/home/acct-mashiyang1991/dhhhh020606/yangchen/knockoff/split.slurm)
    echo "  --> batch jobid ${jobid}"

    # 等待这批全部跑完（RUN/PEND 为 0）
    while true; do
        # 统计真正在跑或排队的 cut_plink 任务数
        nrun=$(squeue -u $USER -n cut_plink -h -t running,pending | wc -l)
        [ "$nrun" -eq 0 ] && break
        echo "[$(date)] Still $nrun tasks running/pending, sleep 60s ..."
        sleep 60
    done
done

echo "[$(date)] All 2495 regions finished."