#!/bin/bash
#SBATCH --job-name=calculate_prs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/prscs/%x%j.log


~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score  /gpfs/commons/home/tlin/output/wightman/prscs/agg_beta.txt  2 4 6 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
--out /gpfs/commons/home/tlin/output/wightman/prscs/prs.


