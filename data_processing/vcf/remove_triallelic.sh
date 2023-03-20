for i in {1..22}
do
## create job
echo create job for chr $i
echo \#!/bin/sh >  remove_triallelic/plink_chr${i}.sh 
echo \#SBATCH --job-name=plink_vcf_chr${i} >>  remove_triallelic/plink_chr${i}.sh
echo \#SBATCH --mail-type=FAIL,END >>  remove_triallelic/plink_chr${i}.sh
echo \#SBATCH --mail-user=tlin@nygenome.org  >>  remove_triallelic/plink_chr${i}.sh
echo \#SBATCH --mem=50G >>  remove_triallelic/plink_chr${i}.sh
echo \#SBATCH --time=10:00:00 >>  remove_triallelic/plink_chr${i}.sh
echo \#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/remove_triallelic/%x.log >>  remove_triallelic/plink_chr${i}.sh

echo >>  remove_triallelic/plink_chr${i}.sh
echo cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink >>  remove_triallelic/plink_chr${i}.sh
echo ~/plink --bfile ADSP_annotated_chr$i --snps-only --exclude ADSP_annotated_whole_genome-merge.missnp --make-bed --out remove_triallelic/chr$i.uni >>  remove_triallelic/plink_chr${i}.sh
sbatch remove_triallelic/plink_chr${i}.sh
done

