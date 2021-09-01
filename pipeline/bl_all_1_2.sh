#!/bin/bash
#SBATCH --job-name=bl_all_1_2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=120G
#SBATCH --time=05:00:00
#SBATCH --output=output/bl_roadmap_deepsea_brain_atac/%x%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
roadmap='/gpfs/commons/home/tlin/polyfun/data/roadmap/roadmap_all_anno.'
brain_atac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr'
deepsea='/gpfs/commons/home/tlin/polyfun/data/deepsea/deepsea_all_anno.'

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix output/bl_roadmap_deepsea_brain_atac/bl_roadmap_deepsea_brain_atac \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$roadmap,$deepsea,$brain_atac \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing


