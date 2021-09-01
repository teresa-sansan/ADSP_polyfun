#!/bin/bash
#SBATCH --job-name=polyfun_1-3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=output/bl_roadmap_deepsea/%x_%j.log


cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

python polyfun.py \
    --compute-ldscores \
    --output-prefix ${output} \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
    --chr ${chr}


