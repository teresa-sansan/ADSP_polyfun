#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/count/%x_%j.log
#SBATCH --array=1-${ttl_blksize}%10
#SBATCH --job-name=extract_bigchunk_chr${chr}

dir_blk="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/"
dir_1kg="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/"
dir_carma='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/'

touch ${dir_carma}/count/count_snp.txt
chr=$1
block=1
set_block=$2
#echo set_block \$SLURM_ARRAY_TASK_ID

for i in /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_ADSP36K_4PRScs_Aug/snplist_ldblk/chr${chr}_*
    do   
        if [ $block -le $set_block ]; then
            actual_snp_count=$(wc -l < $dir_carma/chr${chr}_$block.snplist)
            snp_count_ld=$(wc -l < $i)
            echo $chr $block $snp_count_ld $actual_snp_count >> ${dir_carma}/count/count_snp.txt
            if [ $actual_snp_count -ge 10000 ]; then
                echo $actual_snp_count/count
                cp $i $dir_carma/count
                filename=$(basename $i) 
                mv ${dir_carma}/count/$filename ${dir_carma}/count/big_chunk_blk_${set_block}_${filename}
                echo save file: ${dir_carma}/count/big_chunk_blk${block}_${filename}
            fi
        fi       
        ((block++))
    done


