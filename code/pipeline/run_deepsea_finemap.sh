#!/bin/bash


for i in {1..22}
do	
	echo "run chr$i"
	for j in H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9ac 
	do
        	echo "run $j" 
        	sbatch --export=chr=$i,anno=$j bl_deepsea_finemap.sh
	done
done

