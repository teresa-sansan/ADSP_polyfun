#!/bin/bash
#SBATCH --job-name=intersect_polyfun_prscs_17k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/prscs/original/subset_polyfun/%x%j.log


## note
## this script is using the subset from polyfun (pip > 0) from PRSCS

# prscs='/gpfs/commons/home/tlin/output/wightman/prscs/original/agg_beta.txt'
# subset_polyfun='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_4prscs_not0.tsv'
path='/gpfs/commons/home/tlin/output/wightman/prscs/original/subset_polyfun'
# awk 'NR==FNR { values[$1]=1; next } $2 in values' $subset_polyfun $prscs > $path/interset_polyfun.tsv


~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score $path/interset_polyfun.tsv 2 4 6 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
--out $path/prs
