#!/bin/bash
#SBATCH --job-name=merge_ADSP_Integrate_one_plink
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --partition bigmem
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/%x_%j.log

plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'

for chr in {11..20}
do
	cat $plink_path/36k_ADSP_merge_list.txt| grep chr${chr}.chunk| tail -n+2 > $plink_path/merge_list.tmp


	~/plink	--bfile $plink_path/processing_files/rm_multi_alleleic/ADSP.chr${chr}.chunk1 \
		--merge-list $plink_path/merge_list.tmp \
		--make-bed \
		--allow-no-sex \
		--out $plink_path/ADSP_qc_chr$chr
done
