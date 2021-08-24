#!/bin/bash
#SBATCH --job-name=brain_atac_seq_1-3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=output/kunkle_baseline_microglia_atac/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

python polyfun.py \
    --compute-ldscores \
    --output-prefix output/kunkle_baseline_microglia_atac/kunkle_bl_microglia \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
    --chr ${chr} \
    --allow-missing


