#!/bin/bash
#SBATCH --job-name=remove_multiallelic
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files/rm_multi_alleleic/%x_%j.log

plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'

chunk=1

~/plink	--bfile $plink_path/processing_files/ADSP.chr${chr}.chunk${chunk} \
    --exclude $plink_path/duplist.txt \
	--make-bed \
	--allow-no-sex \
	--out $plink_path/processing_files/rm_multi_alleleic/ADSP.chr${chr}.chunk${chunk}
