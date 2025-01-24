
dir='/gpfs/commons/home/tlin/output/CARMA/'

for anno in bl omics 
do
    cd $dir/$anno
    for i in {1..20}
    do
        touch chr${i}_credibleset.tsv
        for file in ${i}_*.txt.gz
        do
        zcat $file |tail -n+2 | awk '($12 > 0) {print $0}' >> chr${i}_credibleset.tsv
        done
    done
done