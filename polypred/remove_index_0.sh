#!/bin/bash
#SBATCH --job-name=bellenguez_index
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --array=1-22%5
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/remove_index0/%x_%j.log


path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'

#chr=$SLURM_ARRAY_TASK_ID
# for chr in {1..22}
# echo run chr $chr
# for anno in baseline omics omics_dl susie
# do  
#     zcat $path/bellenguez_${anno}_chr${chr}.txt.gz | grep -v '1:0' > $path/remove_index0/bellenguez_${anno}_chr${chr}.txt
# done


for anno in baseline omics omics_dl susie
do  
    cp $path/remove_index0/bellenguez_${anno}_chr1.txt  $path/remove_index0/bellenguez_${anno}.txt
    for chr in {2..22}
    do
        cat $path/remove_index0/bellenguez_${anno}_chr${chr}.txt|tail -n +2 >> $path/remove_index0/bellenguez_${anno}.txt
    done
done
