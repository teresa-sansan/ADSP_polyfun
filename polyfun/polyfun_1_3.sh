#!/bin/bash
#SBATCH --job-name=polyfun_1-3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/%x_%j.log


cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo $output
python polyfun.py \
    --compute-ldscores \
    --output-prefix ${output} \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
    --chr ${chr}


