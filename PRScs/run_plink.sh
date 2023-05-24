#!/bin/bash
#SBATCH --job-name=calculate_prs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/%x%j.log

path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/'
#cat $path/_pst_eff_a1_b0.5_phi1e-02_chr* > $path/agg_prscs_beta.txt

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/agg_prscs_beta.txt 2 4 6 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
--out /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prscs


## compare with directly running with plink (without PRSCS)
if false; then
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score  $path/agg_extract_0.3.tsv  2 4 12 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
--out /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prs_plink

fi
