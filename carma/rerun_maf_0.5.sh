#!/bin/bash
MAX_JOBS=30
check_jobs() {
    # Count the number of currently running sbatch jobs for the user
    #while [ "$(squeue -u tlin | grep 'carma_re' | wc -l)" -ge "$MAX_JOBS" ]; do
    while [ "$(squeue -u tlin | grep 'carma' | wc -l)" -ge "$MAX_JOBS" ]; do
        echo "Waiting for jobs to finish..."
        sleep 5m
    done
}

##Read the input file line by line
# while read -r col1 col2 _; do
#    # Check if the file does not exist
#    if [ ! -e "/gpfs/commons/home/tlin/output/CARMA/geno_filt/remove_maf_0.5/${col1}_${col2}.txt.gz" ]; then
#        echo "run $col1 $col2"
#        check_jobs
#        #python get_nan_snp.py $chr $blk
#        #sbatch --mem=20G  --wrap="python get_nan_snp.py $col1 $col2"
#        sbatch ld_carma_remove_maf0.5.sh "$col1" "$col2" 
#    fi
# done < /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/count/has_nan_failed.txt


while read -r col1 col2 _; do
   if [ ! -e "/gpfs/commons/home/tlin/output/CARMA/geno_filt/${col1}_${col2}.txt.gz" ]; then
        ## check again it's failed
        echo "run $col1 $col2"
        check_jobs
        sbatch run_carma.sh "$col1" "$col2" 
   fi
#done < /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/count/no_nan_failed.txt
