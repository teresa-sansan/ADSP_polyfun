## tested the PRS performance if we change the use effect size * PIP as new effect size 
## because summary stats don't provide PIP, we use the result we got from polyfun_pred
 

cd /gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10

save='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/new_beta/kunkle_max_snp_10.aggregate.csv'
echo $(zcat chr1.aggregrate.all.txt.gz | head -1) NEW_BETA | tr " " "\t" > $save

for i in {1..22}
do
echo 'calculate new beta (effect size * PIP) for chr' ${i}.... 
zcat chr${i}.aggregrate.all.txt.gz | tail -n+2| awk -v OFS='\t' '{print $0 "\t" $10*$11}'>> $save ## column 10 and 11 are PIP and BETA_MEAN
done

echo Start compressing the file...
gzip $save

echo Done!
