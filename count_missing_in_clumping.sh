cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping

for i in {1..22}
do
echo chr $i | tee -a genotyping_rate.txt
missing=$(cat clump_chr${i}.log| grep "is missing from the main dataset, and is a top"|wc -l)
cat clump_chr${i}.log | grep 'Total genotyping rate is' | tee -a genotyping_rate.txt
cat clump_chr${i}.log | grep 'variants and' | tee -a genotyping_rate.txt
#echo 'number of missing SNP =' $missing  | tee -a genotyping_rate.txt
echo ' ' | tee -a genotyping_rate.txt

done
