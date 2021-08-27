#!/bin/bash
#SBATCH --job-name=baseline1-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=UKBB/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

python polyfun.py \
    --compute-h2-bins \
    --output-prefix UKBB/output/kunkle_UKBBbaseline \
    --sumstats data/kunkle_2019_teresa/AD_Kunkle_etal_Stage1.parquet \
    --w-ld-chr /gfps/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2UKB/weights.UKB.


