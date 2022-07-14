cd ~/data/biallelic/check

for chr in {8..22}
do
echo start chr $chr ...
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_UKBB
cut -f 2 ADSP_UKBB_only_${chr}.bim | sort | uniq -d > ${chr}_dups
~/plink --bfile ADSP_UKBB_only_${chr} --exclude ${chr}_dups --make-bed --out ${chr}_filt
done
