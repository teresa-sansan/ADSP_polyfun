#!/bin/bash
#SBATCH --job-name=eqtl_susie
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=30:00:00
#SBATCH --array=601-3000%10

##9047

#/gpfs/commons/home/tlin/script/coloc/snp_of_interest.tsv
line=$SLURM_ARRAY_TASK_ID
MAX_JOBS=200
check_jobs() {
    # Count the number of currently running sbatch jobs for the user
    while [ "$(squeue -u tlin | grep 'ENSG0000' | wc -l)" -ge "$MAX_JOBS" ]; do
        echo "Waiting for jobs to finish..."
        sleep 20m
    done
}
# line=1
start=$((line +1))
head -n $start /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/genelist | tail -1 | while IFS=$'\t' read -r gene
do
    check_jobs 
    sbatch --job-name=$gene eqtl_finemap.sh $gene
done

