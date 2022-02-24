#!/bin/bash
#SBATCH --job-name=bellenguez_clump
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base/%x_%j.log

##bellenguçez_QCed
sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_qc.gz'

##kunkle_QCed
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'
#--out /gpfs/commons/home/tlin/output/cT/kunkle/qc/kunkle_clump_chr${chr}

##bellenguez_not_qced
#sumstats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv.gz'

## Pvalue, SNP

echo start chr $chr
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats \
--clump-snp-field SNP \
--clump-field PVALUE \
--out /gpfs/commons/home/tlin/output/cT/bellenguez/qc_on_base/bellenguez_clump_chr${chr}




##before qc
#--bfile ~/data/biallelic/check/${chr}_filt \
#--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/clump_chr${chr}
#--clump /gpfs/commons/home/tlin/data/kunkle_etal_Stage1_info.tsv \



##qc
#--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr${chr} \

