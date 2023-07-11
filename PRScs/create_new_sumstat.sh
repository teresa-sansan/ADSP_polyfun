## README
## This script is using the output from polyfun then run prscs!

#path=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5
path=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/finemap/max_snp_5
#path=/gpfs/commons/home/tlin/output/kunkle/new_anno/all_anno/finemap/max_snp_5
#echo -n "#org snp:  "
#cat $path/aggregate.all.txt | wc -l

# echo "check in repo" $path
# file_name="agg_extract_pip_not0.tsv"
# echo "extracting SNPs that PIP > 0..."
# awk '{if ($11 > 0)  print $0}' $path/aggregate.all.txt > $path/$file_name
# echo -n "#SNPs with PIP > 0: "
# cat $path/$file_name | wc -l
# #cat $path/agg_extract_0.3.tsv cut -f 2,4,5,12,10 > $path/agg_4prscs.tsv ## the order need to be SNP, A1, A2, BETA and P

# echo "creating the right format for PRSCS..."
# awk '{print $2, $4, $5, $12, $10}' $path/$file_name > $path/agg_4prscs.tsv 
# sed -i '1s/BETA_MEAN/BETA/' $path/agg_4prscs.tsv 

# echo Done! 


## this part is to create new sumstat from org sumstat to meet PRSCS's requirement.

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed'
file_name='Bellenguez_et_al_2021_hg37_new_sep20_qc_nodup.tsv'


head -1 $path/$file_name
echo
echo 'the order need to be SNP, A1, A2, BETA and P'
head -1 $path/$file_name| awk '{print $1, $4, $5, $11, $6}'
awk '{print $1, $4, $5, $11, $6}' $path/$file_name > $path/bellenguez_4prscs.tsv 

echo Done!