path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/correct_chr_vcf/filt'
for i in {4..22}
do
echo start chr${i}
#sbatch --mem=90G -c 5 --wrap="bgzip -c ${path}/ADSP_annotated_chr${i}.genofilt.vcf > ${path}/gzip/ADSP_annotated_chr${i}.genofilt.vcf.gz"
sbatch --mem=90G -c 5 --wrap="gzip ${path}/ADSP_annotated_chr${i}.genofilt.vcf > ${path}/gzip/ADSP_annotated_chr${i}.genofilt.vcf.gz"
done