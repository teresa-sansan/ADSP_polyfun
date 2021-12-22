#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1> /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/max_10.log 2>&1

for chr in {1..22}
do
echo start chr $chr
~/plink --bfile ~/data/biallelic/${chr}_filt \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/chr${chr}.aggregrate.all.txt  \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/clumping/clump_max10_chr${chr}
done
