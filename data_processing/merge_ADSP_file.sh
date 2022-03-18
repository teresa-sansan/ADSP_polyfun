#!/bin/bash
#SBATCH --job-name=merge_ADSP_plinkfile
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/%x.log

#~/plink --allow-no-sex \
#	--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/1_filt \
#	--merge-list /gpfs/commons/home/tlin/data/ADSP_merge_list.txt \
#	--make-bed \
#	--out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/filt/ADSP_merged

~/plink --allow-no-sex \
        --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_chr1 \
        --merge-list /gpfs/commons/home/tlin/data/ADSP_qc_merge_list.txt \
        --make-bed \
        --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc/ADSP_qc_merged

