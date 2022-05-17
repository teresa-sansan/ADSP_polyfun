#!/bin/sh
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position

for i in {1..22}
do
  echo start chr $i
  cat ADSP_annotated_hg37_chr${i}.vcf.tsv | awk '$1=="chr'$i'" {print $0}' > correct_vcf/ADSP_annotated_hg37_chr${i}_correct.vcf.tsv 
done



