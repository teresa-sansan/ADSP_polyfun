#!/bin/bash
#SBATCH --job-name=bl_all_1_2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=120G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle_all/%x%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'
output='/gpfs/commons/home/tlin/output/kunkle_all/all_anno_'

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun

##1-2
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

echo finish polyfun1_2

##1-3
python polyfun.py \
    --compute-ldscores \
    --output-prefix $output \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld

echo finish polyfun1_3

##1-4
python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
    	--sumstats $summary_stats \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing

echo finish polyfun1_4
