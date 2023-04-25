#!/bin/bash

#for anno in bl_dl_annotations bl_brain_atac bl
#for anno in enformer no_ml all_except_enformer update_all+enformer
for anno in old_ml
do
	echo $anno
	for i in 1 5 10
	do
		sbatch --export=max_snp=$i,anno=$anno polypred.sh
		#sbatch --export=max_snp=$i polypred.sh
	done
done

    