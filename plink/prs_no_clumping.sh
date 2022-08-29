#!/bin/bash
#SBATCH --job-name=prs_noclumping
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink/%x_%j.log



## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

for sumstat in bellenguez wightman
do 
  echo running $sumstat
  ##agg_name='aggregate.all.txt' ### this is for kunkle
  agg_name='agg_fixed_converge.tsv'
  for chr in {1..22}
  do
   ~/plink \
   --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_${chr} \
   --score /gpfs/commons/home/tlin/output/$sumstat/${sumstat}_fixed_0224/finemap/max_snp_10/${agg_name} 2 4 11 header \
   --q-score-range range_list.txt /gpfs/commons/home/tlin/output/$sumstat/${sumstat}_fixed_0224/finemap/max_snp_10/aggregate.pvalue  \
   --extract /gpfs/commons/home/tlin/output/cT/new_plink/${sumstat}/fixed_0224/qc/chr${chr}.valid.snp \
   --out /gpfs/commons/home/tlin/output/cT/new_plink/$sumstat/fixed_0224/polyfun_beta_no_clump/polyfun_max10_pT_${chr}
   done
done
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/APOE_2SNPs_qc.tsv 3 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_APOE_qc.tsv 1 4 6 \
#--out /gpfs/commons/home/tlin/output/kunkle/kunkle_qc/2SNP_qc.prs
#--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/kunkle/APOE_qc/kunkle_APOE_qc_chr${chr}
