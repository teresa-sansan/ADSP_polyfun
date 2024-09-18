cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA
for chr in {15..22}
do   
    n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
    echo start checking chr${chr}, with total of $n_ld ld blocks....
    touch count/nan_size.txt
    for i in $(seq 1 "$n_ld")
    do
        nan=$(cat chr${chr}_${i}.ld|grep -o nan |wc -l)
        output=/gpfs/commons/home/tlin/output/CARMA/${chr}_${i}.txt.gz
        if [ -e $output ]; then
            result=$(ls -ll $output | awk '{print $6, $7}')
        else {
           result="failed"
        }
        fi
        echo $chr $i $nan $result | tee -a count/nan_size.txt
    done
done