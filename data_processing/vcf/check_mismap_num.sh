
## check position

hg38='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered/only_position/ADSP_annotated_hg38_chr'
hg37='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/ADSP_annotated_hg37_chr'
correct_vcf='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/correct_vcf/ADSP_annotated_hg37_chr'

for i in {1..22}
do

ori=$(cat ${hg38}$i.vcf.tsv|wc -l)
liftover=$(cat ${hg37}$i.vcf.tsv| wc -l)
correct_map=$(cat ${correct_vcf}${i}_correct.vcf.tsv | wc -l)

percent1=$(awk 'BEGIN{printf "%.2f%\n",('$liftover'/'$ori')*100}')
percent2=$(awk 'BEGIN{printf "%.2f%\n",('$correct_map'/'$liftover')*100}')

echo chr $i in hg38 has $ori SNPs, $liftover " ("$percent1")"  left after liftover.
echo Where $correct_map "("$percent2 ")" SNPs are correctly mapped.  
echo

done
