#!/bin/bash
#SBATCH --job-name=merge_ADSP_Integrate_one_plink
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=200:00:00
#SBATCH --partition bigmem
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/%x_%j.log




############

## setting missing variant id because 
## Warning: Multiple chromosomes seen for variant '.'.
## 	--set-missing-var-ids "." \  ; try fixing the error with this tag but still wasn't solved.
## 	tried 	--set-missing-var-ids @:#"." \ ; still, doesn't work
## 	--flip /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated-merge.missnp \
#	--exclude-snp /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/exclude_SNP.txt \
#	--exclude /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/exclude_SNP.txt \
  

~/plink	--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_chr1 \
	--merge-list /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_merge_list.txt \
	--make-bed \
	--allow-no-sex \
	--out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged


# ~/plink	--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated_chr1  \
# 	--merge-list /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_merge_list.txt \
#    --exclude-snp /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated_whole_genome-merge.missnp \
# 	--make-bed \
# 	--allow-no-sex \
# 	--out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated_merged_v2
