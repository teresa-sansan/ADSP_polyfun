#!/bin/bash
#SBATCH --job-name=wightman
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=180G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/%x%j.log


bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman.munged.parquet'
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo

for i in bl bl_dl_annotations bl_brain_atac bl
do
python polyfun.py \
        --compute-h2-bins \
        --output-prefix /gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/$i/$i \
        --sumstats $summary_stats \
        --w-ld-chr $bl/weights.UKB. \
        --allow-missing

done