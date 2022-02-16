#!/bin/bash
#SBATCH --job-name=polyfun1-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/%x_%j.log

cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/bellenguez'
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename_qc.gz'

python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
    	--sumstats $sumstat \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing


