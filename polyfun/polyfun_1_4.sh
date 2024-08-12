#!/bin/bash
#SBATCH --job-name=polyfun1-4
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/bl/bl/%x_%j.log

cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
output='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/bl/bl'
#sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/kunkle_et_al_2021_hg37_ldsc.munged.parquet'
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sep22_new_bellenguez_et_al_2021_hg37_ldsc.munged.parquet'

python polyfun_assertion_error.py \
 	--compute-h2-bins \
    --output-prefix $output \
    --sumstats $sumstat \
    --w-ld-chr $bl/weights.UKB. \
    --allow-missing


