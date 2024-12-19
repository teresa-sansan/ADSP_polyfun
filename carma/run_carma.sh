#!/bin/sh
#SBATCH --job-name=carma_bl
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=60:00:00
#SBATCH --array=1-2%2
#SBATCH --output=/gpfs/commons/home/tlin/output/CARMA/bl/%x_%j.log

module purge   
module load R/4.2.2    
module load gcc/11.2.0
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so
# chr=$1
# ld=$2
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt
# run by ld
# #ld=$SLURM_ARRAY_TASK_ID
# #size_ld=$(wc -l < chr${chr}_${ld}.bim )   
# Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld
# # if [ "$size_ld" -lt 10000 ]; then
# #      Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld
# # else
# #      echo "skip this block because it has > 10000 SNPs (size: $size_ld)"
# # fi

MAX_JOBS=150
check_jobs() {
# Count the number of currently running sbatch jobs for the user
while [ "$(squeue -u tlin | grep '_ld' | wc -l)" -ge "$MAX_JOBS" ]; do
    echo "Waiting for jobs to finish..."
    sleep 5m
done
}
##----------------------------------------------------
# run by chr 
# only take ~30 min to run per ld chunk, can run by chr
chr=$SLURM_ARRAY_TASK_ID
n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)

for ld in $(seq 11 "$n_ld")
do
    size_ld=$(wc -l < chr${chr}_${ld}.bim )
    if [ "$size_ld" -lt 10000 ] && [ ! -e /gpfs/commons/home/tlin/output/CARMA/bl/${chr}_${ld}.txt.gz ]; then
        check_jobs 
        sbatch --job-name=chr${chr}_ld${ld} --output=/gpfs/commons/home/tlin/output/CARMA/bl/%x_%j.log --mem=15G --wrap="Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld"   
        #Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld
        #Rscript /gpfs/commons/home/tlin/script/carma/test_carma.r $chr $ld
    else
        echo "skip this block because it has > 10000 SNPs (size: $size_ld) or existed"
    fi
done
