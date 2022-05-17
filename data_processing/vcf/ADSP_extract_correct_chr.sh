#!/bin/bash
#SBATCH --job-name=check_correct_chr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct/%x%j.log

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37

for i in {10..21}
do
echo start chr$i
#cat ADSP_annotated_chr${i}.vcf.gz| grep -v '##'| awk -v chr_num="chr$i" '{if($1==chr_num) {print$0}}' > correct/ADSP_chr${i}.vcf
sbatch --mem=15G -c 10 --wrap="zcat ADSP_annotated_chr${i}.vcf.gz | bgzip -c > ADSP_annotated_chr${i}.new.vcf.gz  && tabix ADSP_annotated_chr${i}.new.vcf.gz --n_jobs=10"
done


#for i in {4..9}
#do
#echo start chr$i
#zcat ADSP_annotated_chr${i}.vcf.gz| grep -v '##'| awk -v chr_num="chr$i" '{if($1==chr_num) {print$0}}' > correct/ADSP_chr${i}.vcf
#done
