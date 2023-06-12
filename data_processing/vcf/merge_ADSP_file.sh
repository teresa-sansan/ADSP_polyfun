#!/bin/bash
#SBATCH --job-name=merge_ADSP_Integrate_one_plink
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=200:00:00
#SBATCH --partition bigmem
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc/%x_%j.log




############

## setting missing variant id because 
## Warning: Multiple chromosomes seen for variant '.'.
## 	--set-missing-var-ids "." \  ; try fixing the error with this tag but still wasn't solved.
## 	tried 	--set-missing-var-ids @:#"." \ ; still, doesn't work
## 	--flip /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_annotated-merge.missnp \
#	--exclude-snp /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/exclude_SNP.txt \
#	--exclude /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/exclude_SNP.txt \
  

plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc'

~/plink	--bfile $plink_path/ADSP.chr1.chunk1 \
	--merge-list $plink_path/36k_ADSP_merge_list.txt \
	--make-bed \
	--flip $plink_path/ADSP_annotated_merged-merge.missnp \
	--exclude /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/duplist.txt \
	--exclude-snp /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/duplist.txt \
	--allow-no-sex \
	--out $plink_path/ADSP_annotated_merged


