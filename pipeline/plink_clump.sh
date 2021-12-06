#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1> /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/clumping/max_10.log 2>&1

for chr in {2..3}
do
echo start chr $chr
~/plink --bfile ~/data/biallelic/${chr}_filt \
--clump-r2 0.1  \
--clump-kb 250  \
--clump /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/chr${chr}.aggregrate.all.txt  \
--clump-snp-field SNP \
--clump-field P \
--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/clumping/clump_max10_chr${chr}
done
