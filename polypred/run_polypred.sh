#!/bin/bash

#for anno in bl_dl_annotations bl_brain_atac bl
#for anno in susie
#do
	for i in 1 5 
	do
		#sbatch --export=max_snp=$i,anno=$anno polypred.sh
		sbatch --export=max_snp=$i polypred.sh
	done
#done

