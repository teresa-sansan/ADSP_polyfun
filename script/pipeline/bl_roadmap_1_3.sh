#!/bin/bash
#SBATCH --job-name=bl_roadmap_1_3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=output/bl_roadmap/specific_col/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

python polyfun.py \
    --compute-ldscores \
    --output-prefix output/bl_roadmap/bl_roadmap \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
    --chr ${chr} \
    --allow-missing


