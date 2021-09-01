#!/bin/bash

#outputformat=output/jansenetal/bl_jansen


#outputformat=output/kunkle_baseline_microglia_atac/kunkle_bl_microglia

outputformat=output/bl_deepsea/bl_deepsea

cd ~/polyfun

for i in {1..22}
do
sbatch --export=chr=$i,output=$outputformat polyfun_1_3.sh
done

