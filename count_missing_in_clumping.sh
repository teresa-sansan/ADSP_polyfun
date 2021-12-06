# cd ??
touch missing_SNP.txt

for i in {1..22}
do
echo chr $i | tee -a missing_SNP.txt
missing=$(cat clump_max10_chr${i}.log| grep "is missing from the main dataset, and is a top"|wc -l)
cat clump_max10_chr${i}.log | grep 'Total genotyping rate is' | tee -a missing_SNP.txt
cat clump_max10_chr${i}.log | grep 'variants and' | tee -a missing_SNP.txt
echo 'number of missing SNP =' $missing  | tee -a missing_SNP.txt
echo ' ' | tee -a missing_SNP.txt

done
