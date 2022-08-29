#!/bin/bash
#SBATCH --job-name=new_wightman_clump_genomewide
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/new_plink/wightman/fixed_beta/%x_%j.log


## no qc
plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

if false; then
for ADSP in ADSP
do
if [[ $sumstats != "wightman" ]]; then
echo "not wightman"

~/plink \
--bfile $plink_path/${ADSP}_${chr}\
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink/$sumstats/$ADSP/$sumstats

else

echo wightman
~/plink \
--bfile $plink_path/$ADSP/${ADSP}_all_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink/$sumstats/fixed_beta/$ADSP/$sumstats_${ADSP}_${chr}


fi
done
fi


##qc

if true; then
for ADSP in ADSP_qc_all 
do
if [[ $sumstats != "wightman" ]]; then
echo "not wightman"

~/plink \
--bfile $plink_path/$ADSP/${ADSP}_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink/$sumstats/$ADSP/$sumstats/${ADSP}_${chr} 

else

echo wightman
~/plink \
--bfile $plink_path/$ADSP/${ADSP}_${chr} \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump $sumstats_path \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/cT/new_plink/$sumstats/fixed_beta/$ADSP/$sumstats_${ADSP}_${chr}


fi
done

fi
