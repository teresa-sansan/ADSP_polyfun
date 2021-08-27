#!/bin/bash

cd ~/polyfun

for i in {1..22}
do
sbatch --export=chr=$i baseline_microglia_atac_1_3.sh
done

