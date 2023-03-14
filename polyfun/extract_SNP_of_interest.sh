path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/'
path_converge='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/'
p_value=10
if false; then
for max_snp in max_snp_1 max_snp_5
do
    cd $path/$max_snp
    pwd
    zcat chr11.aggregate.all.txt.gz| head -n 1 > agg_extract_1e-3.tsv
    echo "start extracting SNP with p < 1e-3 "

    ### kunke  p value is in 9th column, while bellenguez and wightman are in 10th
    echo 'start extracting pvalue < 0.001'
    for i in {1..22}
    do  
        echo Start chr ${i} 
        zcat chr${i}.aggregate.all.txt.gz|tail -n+2|awk -v p_value_col=$p_value '{if($p_value_col <= 0.001) print$0}' >> agg_extract_1e-3.tsv
    done
    echo finished, total line = $( cat agg_extract_1e-3.tsv| wc -l )
done
fi

if true;then
for anno in bl susie all_anno
do
cd $path_converge/$anno/finemap/max_snp_10
zcat agg_fixed_converge.tsv.gz| tail -n+2|awk -v p_value_col=$p_value '{if($p_value_col <= 0.001) print$0}' >> agg_extract_1e-3_fix_converge.tsv
done
fi
