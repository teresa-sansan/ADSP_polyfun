#!/bin/sh 
#SBATCH --job-name=liftover_variant_counts
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/%x_%j.log


## check position
output='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/map_result_mismap_chr${i}.txt'
hg38='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered/only_position/ADSP_annotated_hg38_chr'
hg37='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/ADSP_annotated_hg37_chr'
wrong_vcf='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/incorrect_vcf/new_ADSP_annotated_hg37_chr'
correct_vcf='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct/ADSP_annotated_chr'


echo checking chr ${i} ....
ori=$(cat ${hg38}$i.vcf.tsv|wc -l)
liftover=$(cat ${hg37}$i.vcf.tsv| wc -l)
correct_map=$(cat ${correct_vcf}${i}.correct.vcf.gz | cut -f 1|wc -l)
wrong=$(cat ${wrong_vcf}${i}_incorrect.vcf.tsv | wc -l )

percent1=$(awk 'BEGIN{printf "%.2f%\n",('$liftover'/'$ori')*100}')
percent2=$(awk 'BEGIN{printf "%.2f%\n",('$correct_map'/'$liftover')*100}')
percent3=$(awk 'BEGIN{printf "%.2f%\n",('$wrong'/'$liftover')*100}')

echo chr $i in hg38 has $ori SNPs, $liftover " ("$percent1")"  left after liftover. |tee -a  $output
echo Where $correct_map "("$percent2 ")" SNPs are correctly mapped in the same chromosome. | tee -a $output
echo There are $wrong "("$percent3 ")" SNPs incorrectly mapped in the chr $i. | tee -a $output
echo |tee -a $output


