#!/bin/bash

#for anno in bl_dl_annotations bl_brain_atac bl
#for anno in susie
#do
	for i in 1 5 10
	do
		#sbatch --export=max_snp=$i,anno=$anno polypred.s sh
		sbatch --export=max_snp=$i polypred.sh
	done
#done

    