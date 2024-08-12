#!/bin/bash
#SBATCH --job-name=sbayesc_score
#SBATCH --partition=pe2
#SBATCH --nodes=1           # minimum number of nodes to be allocated
#SBATCH --ntasks=1          # number of tasks
#SBATCH --cpus-per-task=8   # number of cores on the CPU for the task
#SBATCH --mem=50G
#SBATCH --time=19:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-22%13
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/%x_%j.log

chr=$SLURM_ARRAY_TASK_ID

output='/gpfs/commons/home/tlin/output/sbayesRC/bellenguez_whole_genome/bellenguez_tune'
# ~/plink \
#     --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr} \
#     --score ${output}.sbrc.txt 1 2 3 header \
#     --out ${output}_chr${chr}



awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' ${output}_*.profile > ${output}.prs
echo wrote  ${output}.prs
