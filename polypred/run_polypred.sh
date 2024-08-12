#!/bin/bash

#for anno in bl_dl_annotations bl_brain_atac bl
#for anno in enformer no_ml all_except_enformer update_all+enformer
#for anno in no_ml bl all_enformer all_anno 
#for anno in all_except_enformer bl no_ml_new
for anno in only_ml
do
	echo $anno 
	for i in  5
	do
		sbatch --export=max_snp=$i,anno=$anno polypred.sh
		#sbatch --export=max_snp=$i polypred.sh
	done
done

    