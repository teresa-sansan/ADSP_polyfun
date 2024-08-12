#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'
path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'

for anno in bl
do
    cd $path/$anno/finemap
    echo start $anno
    for i in 1 3 5 10
    do
        echo start in max_snp_$i
        cat max_snp_${i}/aggregate.all.txt| head -1 | cut -f 10
        cat max_snp_${i}/aggregate.all.txt| head -1 > max_snp_${i}/agg_extract_1e-3.tsv
        cat max_snp_${i}/aggregate.all.txt |awk '{if($10 <= 0.001) print$0}' >> max_snp_${i}/agg_extract_1e-3.tsv

    done
    echo
done