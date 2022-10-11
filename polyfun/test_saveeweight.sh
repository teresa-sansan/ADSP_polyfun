#!/bin/bash
#SBATCH --job-name=bl_all_1_2
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=120G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/test_saveweight.log

python /gpfs/commons/home/tlin/polyfun_omer_repo/polyfun_save_weight.py --compute-h2-L2 --output-prefix /gpfs/commons/home/tlin/output/test0909 --sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet --ref-ld-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.,/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K27ac/merged_annotations_ukb/brain_H3K27ac_seq_chr --w-ld-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/weights.UKB. --allow-missing
