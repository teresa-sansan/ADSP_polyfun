#!/bin/sh
#SBATCH --job-name=carma
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=190G
#SBATCH --time=60:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/CARMA/%x_%j.log
#SBATCH --array=10-20%5

module purge   
module load R/4.2.2    
module load gcc/11.2.0
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA
# run by ld
# chr=22
# ld=$SLURM_ARRAY_TASK_ID
# 

# size_ld=$(wc -l < chr${chr}_${ld}.bim )
    
# if [ "$size_ld" -lt 10000 ]; then
#      Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld

# else
#      echo skip this block because it has > 10000 SNPs (size: $size_ld)
# fi


##----------------------------------------------------
## run by chr 
## only take ~30 min to run per ld chunk, can run by chr
chr=$SLURM_ARRAY_TASK_ID

n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
for ld in $(seq 1 "$n_ld")
do
    size_ld=$(wc -l < chr${chr}_${ld}.bim )
    if [ "$size_ld" -lt 10000 ] && [ ! -e /gpfs/commons/home/tlin/output/CARMA/${chr}_${ld}.txt.gz ]; then
        Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $ld
    else
        echo "skip this block because it has > 10000 SNPs (size: $size_ld) or existed"
    fi
done