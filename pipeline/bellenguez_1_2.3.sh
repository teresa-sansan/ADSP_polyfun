#!/bin/bash
#SBATCH --job-name=bellenguez1_2.3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=249G
#SBATCH --time=08:00:00
#SBATCH --output=output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/%x_%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/'
brain_atac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr'
roadmap='/gpfs/commons/home/tlin/polyfun/data/roadmap/roadmap_all_anno.'
deepsea='/gpfs/commons/home/tlin/polyfun/data/deepsea/deepsea_all_anno.'
output='output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac'

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$roadmap,$brain_atac,$deepsea \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing 

python polyfun.py \
  --compute-ldscores \
  --output-prefix $output \
  --ld-ukb
  --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld \
  --allow-missing 
 



