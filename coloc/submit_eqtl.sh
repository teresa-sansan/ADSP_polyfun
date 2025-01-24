#!/bin/bash
#SBATCH --job-name=eqtl_susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --array=1-17%5
#SBATCH --output=/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/%x_%j.log

##9000-9047
##4000-8999
##1000-4109
line=$SLURM_ARRAY_TASK_ID
MAX_JOBS=200
check_jobs() {
    # Count the number of currently running sbatch jobs for the user
    while [ "$(squeue -u tlin | grep 'ENSG0000' | wc -l)" -ge "$MAX_JOBS" ]; do
        echo "Waiting for jobs to finish..."
        sleep 5m
    done
}
# line=1
start=$((line +1))
#head -n $start /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/genelist | tail -1 | while IFS=$'\t' read -r gene
#head -n $start /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene3.txt | tail -1 | while IFS=$'\t' read -r gene

head -n $start /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss_rerun_gene2.txt | tail -1 | while IFS=$'\t' read -r gene
do
    check_jobs 
    if [ ! -e /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_*${gene}.rds ];then
        sbatch --job-name=$gene eqtl_finemap.sh $gene 
    fi
done

