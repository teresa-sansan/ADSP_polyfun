## This is fast. Don't even need to submitted to cluster. 

## kunkle
if false; then
sumstats_path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10'
--score $sumstats_path/aggregate.all.txt 2 4 11 header 
fi

##bellenguez
if true; then
#sumstats_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10'
sumstats_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/susie/max_snp_10'
agg_name='aggregate.all.txt'
sumstat='bellenguez'
#--score $sumstats_path/agg_fixed_converge.tsv 2 4 11 header \
#--extract $path/$sumstat/fixed_0224/polyfun_beta/${chr}.valid.snp \
fi

##wightman
if false; then
sumstats_path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/'
sumstat='wightman'
--score $sumstats_path/agg_fixed_converge.tsv 2 4 11 header \

fi

plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink'

echo "creating non-zipped sumstat..."
#zcat $sumstats_path/${agg_name} > $sumstats_path/${agg_name}.tsv
cat  $sumstats_path/${agg_name} | cut -f 2,9 > $sumstats_path/aggregate.pvalue
path='/gpfs/commons/home/tlin/output/cT/new_plink/'

echo "start p value thresholding."
for chr in {1..22}
do
#awk 'NR!=1{print $3}' $path/$sumstat/fixed_0224/polyfun_beta/${sumstat}_polyfun_qc_chr${chr}.clumped > $path/$sumstat/fixed_0224/polyfun_beta/${chr}.valid.snp

~/plink \
--bfile $plink_path/ADSP_qc_all/ADSP_qc_all_${chr} \
--score $sumstats_path/${agg_name} 2 4 11 header \
--q-score-range range_list.txt $sumstats_path/aggregate.pvalue  \
--extract /gpfs/commons/home/tlin/output/cT/new_plink/wightman/fixed_0224/qc/chr${chr}.valid.snp \
--out $path/$sumstat/fixed_0224/polyfun_beta/polyfun_max10_susie_pT_${chr} 


done

echo 'done'
