#!/bin/bash 

## kunkle
sumstat='kunkle_et_al_2021_hg37_ldsc.munged.parquet'
output='kunkle'
#sbatch  --export=path_sumstat=$sumstat,path_out=$output calculate_tau_se.sh

## bellenguez
sumstat='Bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
output='bellenguez'
#sbatch  --export=path_sumstat=$sumstat,path_out=$output calculate_tau_se.sh


## wightman
sumstat='Wightman.munged.parquet'
output='wightman'
sbatch  --export=path_sumstat=$sumstat,path_out=$output calculate_tau_se.sh


## jansen
sumstat='jansen_et_al_2021_hg37_ldsc.munged.parquet'
output='jansen'
#sbatch  --export=path_sumstat=$sumstat,path_out=$output calculate_tau_se.sh
