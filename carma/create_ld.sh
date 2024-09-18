#!/bin/bash

chr=$1
ttl_blksize=$2
tmp_script=$(mktemp)

cat <<EOT > $tmp_script
#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/%x_%j.log
#SBATCH --array=1-${ttl_blksize}%10
#SBATCH --job-name=create_ld_carma_chr${chr}

dir_blk="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/"
dir_1kg="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/"

block=1
echo set_block \$SLURM_ARRAY_TASK_ID

for i in /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_ADSP36K_4PRScs_Aug/snplist_ldblk/chr${chr}_*
    do   
        if [ ! -e ldblk/chr${chr}_\${block}.ld ] && ((block == SLURM_ARRAY_TASK_ID )) ; then
            echo run blk \$block
            /gpfs/commons/home/tlin/plink \
            --bfile \${dir_1kg}/ADSP_EUR_chr${chr} \
            --keep-allele-order --extract \$i --r --matrix --out \${dir_blk}/chr${chr}_\${block} --memory 2000 --write-snplist --maf --make-just-bim --freq --snps-only --geno
            echo created ld
        fi       
        ((block++))
    done
EOT

# Submit the temporary script
sbatch $tmp_script

# Optionally remove the temporary script after submission
rm $tmp_script
