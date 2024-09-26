#!/bin/bash
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/remove_maf_0.5/%x_%j.log
#SBATCH --job-name=rerun_maf_0.5

dir_blk="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/remove_maf_0.5"
dir_1kg="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/"

chr=$1
blk=$2

##create ld
echo run chr$chr, blk$blk
cat $dir_blk/../count/chr${chr}_${blk}_nan_snp.tsv | tail -n+2| cut -f 2 > $dir_blk/../count/chr${chr}_${blk}_remove_snp.txt

/gpfs/commons/home/tlin/plink \
--bfile ${dir_1kg}/ADSP_EUR_chr${chr} \
--keep-allele-order --extract $dir_blk/../chr${chr}_${blk}.snplist --exclude $dir_blk/../count/chr${chr}_${blk}_remove_snp.txt \
--r --matrix --out ${dir_blk}/chr${chr}_${blk} --memory 2000 \
--write-snplist --maf --make-just-bim --freq --snps-only --geno

## run carma
module purge   
module load R/4.2.2    
module load gcc/11.2.0
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt
Rscript /gpfs/commons/home/tlin/script/carma/carma.r $chr $blk True