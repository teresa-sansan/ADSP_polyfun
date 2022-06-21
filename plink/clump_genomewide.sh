#!/bin/bash
#SBATCH --job-name=clump_genomewide
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/genomewide_plink/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

if true; then
for ADSP in ADSP_all ADSP_UKBB_only 
do
if [[ $sumstats != "wightman" ]]; then
echo "not wightman"

~/plink \
--bfile $plink_path/$ADSP \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/$plink/$sumstats

else

echo wightman
~/plink \
--bfile $plink_path/$ADSP \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field Pvalue \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/$plink/$sumstats


fi
done
fi


##qc

if false; then
for ADSP in ADSP_all_qc ADSP_UKBB_qc
do
if [[ $sumstats != "wightman" ]]; then
echo "not wightman"

~/plink \
--bfile $plink_path/$ADSP \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/$plink/$sumstats

else

echo wightman
~/plink \
--bfile $plink_path/$ADSP \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/$sumstats/$plink/$sumstats


fi
done

fi
