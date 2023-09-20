#!/bin/bash
#SBATCH --job-name=prs_17k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSCS/17k/bellenguez/%x%j.log

#path='/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer'
#path='/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/beta_sumstat'
path='/gpfs/commons/home/tlin/output/prs/PRSCS/17k/bellenguez'
cat $path/_pst_eff_a1_b0.5_phi1e-02_chr* > $path/agg_prscs_beta.txt

bellenguez_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc.pvalue'
wightman_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp'

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score $path/agg_prscs_beta.txt 2 4 6 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt $bellenguez_snp \
--out $path/prs
