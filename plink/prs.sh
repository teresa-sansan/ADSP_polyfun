#!/bin/bash
#SBATCH --job-name=bellenguez_adsp_reference
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --array=1-22%12
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/%x_%j.log

## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

polyfun_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
sumstat='bellenguez'

chr=$SLURM_ARRAY_TASK_ID

for anno in susie 
do 
  echo running $anno
    zcat $polyfun_path/${sumstat}_${anno}_chr${chr}.txt.gz > $output/bellenguez_adsp_reference/${anno}_chr${chr}.txt
    ~/plink \
    --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink/ADSP_chr${chr}_filt \
    --score $output/bellenguez_adsp_reference/${anno}_chr${chr}.txt 1 4 11 header \
    --out $output/bellenguez_adsp_reference/thres/${anno}_chr${chr}

    rm $output/bellenguez_adsp_reference/${anno}_chr${chr}.txt
  
done

#  susie: 1,4,11
## others: 2,4,12 baseline omics omics_dl