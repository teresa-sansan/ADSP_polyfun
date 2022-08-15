#!/bin/bash
#SBATCH --job-name=pT_kunkle_polyfun
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink/kunkle/fixed_0224/polyfun_beta/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

bfile='ADSP_qc_all'

path='/gpfs/commons/home/tlin/output/cT/new_plink/kunkle/fixed_0224/polyfun_beta'
awk 'NR!=1{print $3}' $path/kunkle_polyfun_qc_chr${chr}.clumped > $path/${chr}.valid.snp

~/plink \
--bfile $plink_path/${bfile}_${chr} \
--score  $sumstats_path 2 4 $BETA_MEAN header \
--q-score-range range_list.txt $path/aggregate.pvalue  \
--extract $path/${chr}.valid.snp \
--out $path/${bfile}_pT_${chr} 


