#!/bin/sh
#SBATCH --job-name=liftover_correction
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/vcf_correct/%x_%j.log

#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/


echo start chr $i
#cat ADSP_annotated_hg37_chr${i}.vcf.tsv | awk '$1=="chr'$i'" {print $0}' > correct_vcf/ADSP_annotated_hg37_chr${i}_correct.vcf.tsv 

bcftools view -r chr$i annotated_hg37/ADSP_chr${i}.vcf.gz -O z -o  annotated_hg37_plink_ibd/vcf_correct/ADSP_chr${i}.vcf

