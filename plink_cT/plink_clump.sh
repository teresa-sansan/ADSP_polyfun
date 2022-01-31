#!/bin/bash
#SBATCH --job-name=kunkle_clump
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/cT/kunkle/%x_%j.log



echo start chr $chr

~/plink \
--bfile ~/data/biallelic/check/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/home/tlin/data/kunkle_etal_Stage1_info.tsv \
--clump-snp-field SNP \
--clump-field Pvalue \
--out /gpfs/commons/home/tlin/output/cT/kunkle/clump_chr${chr}

#--clump /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/chr${chr}.aggregrate.all.txt  \
#--clump /gpfs/commons/home/tlin/data/bellenguez_2021_final.tsv.gz \
#--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/clump_chr${chr}
