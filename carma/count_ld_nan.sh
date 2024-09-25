cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt
for chr in {1..22}
do   
    n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
    echo start checking chr${chr}, with total of $n_ld ld blocks....
    touch count/nan_size.txt
    for i in $(seq 1 "$n_ld")
    do
        size_ld=$(wc -l < chr${chr}_${i}.bim)
        nan=$(cat chr${chr}_${i}.ld|grep -o nan |wc -l)
        output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/${chr}_${i}.txt.gz
        if [ -e $output ]; then
            result='success'
        else {
           result="failed"
        }
        fi
        echo $chr $i $size_ld $nan $result | tee -a count/nan_size.txt
    done
done

