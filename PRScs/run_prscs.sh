#!/bin/bash

## specify what to run
bellenguez='bellenguez_4prscs.tsv'
n_bellenguez=487511

wightman='wightman_prscs.tsv'
n_wightman=762971

output_dir='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_ibd_adsp_fixed/'

for chr in 22
do
	#sbatch --export=chr=$chr,n_sumstat=$n_bellenguez,sumstat=$bellenguez,output_dir=$output_dir/bellenguez prscs.sh
	sbatch --export=chr=$chr run_plink_36k.sh
done