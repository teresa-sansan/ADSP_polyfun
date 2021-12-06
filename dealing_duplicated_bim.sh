cd ~/data/biallelic/

for chr in {2..3}
do
echo start chr $chr ...
cut -f 2 /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/ADSP_chr${chr}.bim | sort | uniq -d > ~/data/biallelic/${chr}_dups
~/plink --bfile ADSP_chr${chr} --exclude ${chr}_dups --make-bed --out ${chr}_filt
done
