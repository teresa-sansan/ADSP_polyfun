
path=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5
echo -n "#org snp:  "
#cat $path/aggregate.all.txt | wc -l

echo

echo -n "#SNPs with PIP >= 0.3: "
# awk '{if ($11 >= 0.3)  print $0}' aggregate.all.txt > $path/agg_extract_0.3.tsv
cat $path/agg_extract_0.3.tsv | wc -l
#cat $path/agg_extract_0.3.tsv cut -f 2,4,5,12,10 > $path/agg_4prscs.tsv ## the order need to be SNP, A1, A2, BETA and P

awk '{print $2, $4, $5, $12, $10}' $path/agg_extract_0.3.tsv > $path/agg_4prscs.tsv 
sed -i '1s/BETA_MEAN/BETA/' $path/agg_4prscs.tsv 


echo  -n "which are in the following CHR:"
cat $path/agg_extract_0.3.tsv| tail -n+2 | cut -f 1| uniq| tr '\n' ','

echo
echo







#awk '{if ($11 > 0.3)  print $0}' aggregate.all.txt > agg_extract_0.3.tsv