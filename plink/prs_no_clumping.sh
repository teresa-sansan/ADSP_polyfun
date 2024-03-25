#!/bin/bash
#SBATCH --job-name=prs_bellenguez
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=10:00:00
#SBATCH --array=1-22%12
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/%x_%j.log

## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

polyfun_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fine_mapping/annotations_dl/aggregate_finemap/'
adsp_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/'
sumstat='bellenguez'

chr=$SLURM_ARRAY_TASK_ID
thres=0.25
for anno in baseline omics omics_dl 
do 
  echo running $anno
    zcat $adsp_path/${sumstat}_${anno}_chr${chr}.txt.gz > $output/bellenguez_adsp_reference/thres/${sumstat}_${anno}_chr${chr}.txt
    
    ~/plink \
    --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr} \
    --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/pip_filt/${sumstat}_${anno}_pip_${thres}.snp \
    --score  $output/bellenguez_adsp_reference/thres/${sumstat}_${anno}_chr${chr}.txt 2 4 12 header \
    --out $output/bellenguez_adsp_reference/thres/${anno}_pip${thres}_chr${chr}
    
    rm  $output/bellenguez_adsp_reference/thres/${sumstat}_${anno}_chr${chr}.txt
  
done
