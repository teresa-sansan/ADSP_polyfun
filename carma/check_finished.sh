anno=$1
#path='/gpfs/commons/home/tlin/output/CARMA/$anno/'
cd /gpfs/commons/home/tlin/output/CARMA/$anno/
touch ../${anno}_unfinished3.txt

for chr in {1..22}
do
    echo running on chr$chr
    n_ld=$(ls -1 /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt | grep .bim | grep chr${chr}_ | wc -l)
    for ld in $(seq 1 "$n_ld")
    do 
    if [ ! -e ${chr}_${ld}.txt.gz ]; then
        echo ${chr} ${ld}  >> ../${anno}_unfinished3.txt
    fi
    done
done
