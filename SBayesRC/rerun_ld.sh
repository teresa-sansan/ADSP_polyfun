#!/bin/bash
#SBATCH --job-name=chr22_gctb   
#SBATCH --mem=16G           
#SBATCH --array=2088-2108%10
#SBATCH --time=02:00:00            
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/%x_%j.log.log

source /gpfs/commons/home/tlin/miniconda3/etc/profile.d/conda.sh
conda activate sbayesrc
module load R/4.3.3
chr=22


## this step regenerate the bim (remove the duplicated SNP)
## rs75130031 in chr7 and rs116652232 in chr19
## only need to be ran once

# ~/plink --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr19 \
#         --exclude /gpfs/commons/home/tlin/script/SBayesRC/duplicated.snp \
#         --make-bed --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/no_dup_ADSP_EUR_chr19

blk=${SLURM_ARRAY_TASK_ID}

# gctb --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr${chr} \
# --chr ${chr} --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/snplist/${blk}.snplist --make-full-ldm \
# --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/b${blk}

# cp  /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/b${blk}* /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/LD_sbayesrc/

gctb --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/plink_file/ADSP_EUR_chr22 --chr 22 --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/snplist/${blk}.snplist --make-full-ldm --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/b2088 &> /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc/b{blk}.log
