cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt
echo chr ld_blk ld_size nan result remove_maf_0.5 > count/nan_size_rerun2.txt
for chr in {1..22}
do   
    n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
    echo start checking chr${chr}, with total of $n_ld ld blocks....
    for i in $(seq 1 "$n_ld")
    do
        if [ -e remove_maf_0.5/chr${chr}_${i}.bim ]; then
            size_ld=$(wc -l < remove_maf_0.5/chr${chr}_${i}.bim)
            nan=$(cat remove_maf_0.5/chr${chr}_${i}.ld|grep -o nan |wc -l)
            output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/remove_maf_0.5/${chr}_${i}.txt.gz
            remove_maf='True'

        else
            size_ld=$(wc -l < chr${chr}_${i}.bim)
            nan=$(cat chr${chr}_${i}.ld|grep -o nan |wc -l)
            output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/${chr}_${i}.txt.gz
            remove_maf='False'
        fi

        size_ld=$(wc -l < chr${chr}_${i}.bim)
        nan=$(cat chr${chr}_${i}.ld|grep -o nan |wc -l)
        output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/${chr}_${i}.txt.gz
        if [ -e $output ]; then
            result='success'
        else {
           result="failed"
        }
        fi
        echo $chr $i $size_ld $nan $result $remove_maf | tee -a count/nan_size_rerun2.txt
    done
done

