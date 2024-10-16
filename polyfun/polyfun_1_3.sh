#!/bin/bash
#SBATCH --job-name=polyfun_1-3
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0824/%x_%j.log


cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo $output

python polyfun.py \
    --compute-ldscores \
    --output-prefix ${output} \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
    --chr $chr
