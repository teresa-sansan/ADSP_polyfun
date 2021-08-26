#!/bin/bash
#SBATCH --job-name=bl_brain_atac_1_2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=80G
#SBATCH --time=-05:00:00
#SBATCH --output=output/kunkle_UKBB_brain_atac/%x%j.log

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
brain_atac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb'


python polyfun.py \
  --compute-h2-L2 \
  --output-prefix output/kunkle_UKBB_brain_atac/kunkle_brain_atac \
  --sumstats data/kunkle_2019_teresa/AD_Kunkle_etal_Stage1.parquet \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$brain_atac/brain_atac_seq_chr \
  --w-ld-chr data/data_backup/baseline_backup/weights.UKB. \
  --allow-missing


