#!/bin/sh
## check position
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position
extract_vcf='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/ADSP_annotated_hg37_chr'

for i in {1..4}
do
cat ${extract_vcf}${i}.vcf.tsv | awk '$1 != "chr'$i'"{print $0}' > incorrect_vcf/new_ADSP_annotated_hg37_chr${i}_incorrect.vcf.tsv &
done
