#!/bin/bash
#SBATCH --job-name=bellenguez
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/%x%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'
brain_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K4me3/merged_annotations_ukb/brain_H3K4me3_seq_chr'
brain_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K27ac/merged_annotations_ukb/brain_H3K27ac_seq_chr'
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
#output='/gpfs/commons/home/tlin/output/kunkle_all_2/all_anno'


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo

##1-2
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno,$brain_H3K4me3,$brain_H3K27ac \
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
