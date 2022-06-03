for i in {1..22}
do
## create job
echo create job for chr $i
echo \#!/bin/sh >  plink_file_job/plink_chr${i}.sh 
echo \#SBATCH --job-name=plink_vcf_chr${i} >>  plink_file_job/plink_chr${i}.sh
echo \#SBATCH --mail-type=FAIL,END >>  plink_file_job/plink_chr${i}.sh
echo \#SBATCH --mail-user=tlin@nygenome.org  >>  plink_file_job/plink_chr${i}.sh
echo \#SBATCH --mem=50G >>  plink_file_job/plink_chr${i}.sh
echo \#SBATCH --time=10:00:00 >>  plink_file_job/plink_chr${i}.sh
echo \#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/keep_biallelic/%x.log >>  plink_file_job/plink_chr${i}.sh

echo >>  plink_file_job/plink_chr${i}.sh
echo cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink >>  plink_file_job/plink_chr${i}.sh
echo ~/plink --bfile ADSP_annotated_chr$i --snps-only --exclude ADSP_annotated-merge.missnp --make-bed --out keep_biallelic/chr$i.uni >>  plink_file_job/plink_chr${i}.sh
sbatch plink_file_job/plink_chr${i}.sh
done

