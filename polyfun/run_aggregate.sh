#!/bin/bash
#bl_brain_atac bl_dl_annotations
#for anno in susie all_anno
#for anno in bl_brain_atac bl_dl_annotations bl
#for anno in susie
#do
for i in  1  10
do
        echo "max_snp = $i" 
        #sbatch --export=max_snp=$i aggregate_finemap_result.sh 
        sbatch --export=max_snp=$i,anno=$anno  aggregate_finemap_result.sh 
	#sbatch --export=max_snp=$i aggregate_finemap_result_kunkle.sh
done
#done
