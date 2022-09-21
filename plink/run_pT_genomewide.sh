#!/bin/bash

for i in {1..21}
do
sbatch --export=chr=$i,sumstats='wightman/fixed_beta',sumstats_path='wightman_fixed_beta_nodup.tsv',pvalue='wightman_fixed_beta.pvalue',beta=10 pT.sh
done
