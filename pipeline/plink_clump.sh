#!/bin/bash
#SBATCH --job-name=bellenguez_clump
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/%x_%j.log



echo start chr $chr
~/plink \
--bfile ~/data/biallelic/${chr}_filt \
--clump-p1 1 \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/home/tlin/data/bellenguez_2021_final.tsv.gz \
--clump-snp-field RSID \
--clump-field PVALUE \
--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/clump_chr${chr}


#--clump /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/chr${chr}.aggregrate.all.txt  \
