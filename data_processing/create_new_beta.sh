## Teresa Lin
## Last edit: 07.12.22 
## Motivation: tested the PRS performance if we change the use effect size * PIP as new effect size 
## because summary stats don't provide PIP, we use the result we got from polyfun_pred
 

##kunkle
if false; then
cd /gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10
save='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/new_beta/kunkle_max_snp_10.aggregate.tsv'
echo $(zcat chr1.aggregrate.all.txt.gz | head -1) NEW_BETA | tr " " "\t" > $save
fi

##wightman

#save='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/wightman_agg_fixed_converge_new_beta.tsv'
#file='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/agg_fixed_converge.tsv.gz'
#echo $(zcat $file | head -1) NEW_BETA | tr " " "\t" > $save   ## create column for new file. 

for max_snp in max_snp_1 max_snp_3 max_snp_5 max_snp_7 max_snp_10
do
 echo 'Start processing in max_snp_' $max_snp
 cd /gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/$max_snp
 pwd
 save=$(echo "/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/"${max_snp}"/agg_all_new_beta.tsv")
 echo $(zcat chr1.aggregrate.all.txt.gz | head -1) NEW_BETA | tr " " "\t" > $save


## run this when you are still running chr-separated file. 
 if true; then
 for i in {1..22}
 do
 echo 'calculate new beta (effect size * PIP) for chr' ${i}.... 
 zcat chr${i}.aggregrate.all.txt.gz | tail -n+2| awk -v OFS='\t' '{print $0 "\t" $10*$11}'>> $save ## column 10 and 11 are PIP and BETA_MEAN
 done
 fi


## run this if they are already aggregate into one file.
if false; then 
echo 'Start calculating new beta (effect size * PIP) for' $file
zcat $file | tail -n+2 | awk -v OFS='\t' '{print $0 "\t" $10*$11}'>> $save ## column 10 and 11 are PIP and BETA_MEAN
echo 'Finished, result wrote in ' $save
fi

echo Start compressing the result...
gzip $save

echo Done!
done
