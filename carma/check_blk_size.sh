cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA
for chr in {1..22}
do
   
    n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
    echo start checking chr${chr}, with total of $n_ld ld blocks....
    touch count/blksize.txt

    for i in $(seq 1 "$n_ld")
    do
        size_ld=$(wc -l < chr${chr}_${i}.bim )
        
        echo $chr $size_ld >> count/blksize.txt

    done
done