cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37

for i in {10..22}
do
echo start chr$i
cat ADSP_annotated_chr${i}.vcf.gz| grep -v '##'| awk -v chr_num="chr$i" '{if($1==chr_num) {print$0}}' > correct/ADSP_chr${i}.vcf
done


for i in {4..9}
do
echo start chr$i
zcat ADSP_annotated_chr${i}.vcf.gz| grep -v '##'| awk -v chr_num="chr$i" '{if($1==chr_num) {print$0}}' > correct/ADSP_chr${i}.vcf
done
