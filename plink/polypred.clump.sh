#!/bin/bash
#SBATCH --job-name=polyfun_beta_pt_wightman
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink/wightman/fixed_0224/polyfun_beta/%x_%j.log

##bellenguez_QCed
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_qc.tsv.gz'

##belenguez_munged
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_ldsc.munged.parquet'

##kunkle_QCed
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'
#--out /gpfs/commons/home/tlin/output/cT/kunkle/qc/kunkle_clump_chr${chr}

#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv.gz'
output='/gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224'
pvalue="P"
##wightman
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_qc.tsv.gz'
#output='/gpfs/commons/home/tlin/output/cT/wightman'
## Pvalue, SNP

echo sumstat= $sumstat

echo start chr $chr
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink/bellenguez/fixed_0224/polyfun_beta/bellenguez_polyfun_qc_chr${chr}


