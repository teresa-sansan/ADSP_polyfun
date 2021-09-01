#!/bin/bash
#SBATCH --job-name=bl_roadmap_microglia_1_2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=400G
#SBATCH --time=05:00:00
#SBATCH --output=output/bl_roadmap_microglia/%x_%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
brain_atac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr'

roadmap_DNase='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr'
roadmap_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K27ac/roadmap_H3K27ac_chr'
roadmap_H3K4me1='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me1/roadmap_H3K4me1_chr'
roadmap_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me3/roadmap_H3K4me3_chr'


python polyfun.py \
  --compute-h2-L2 \
  --output-prefix output/bl_roadmap_microglia/bl_roadmap_microglia \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$roadmap_DNase,$roadmap_H3K27ac,$roadmap_H3K4me1,$roadmap_H3K4me3,$brain_atac  \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing 
