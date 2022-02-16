#!/bin/bash
#SBATCH --job-name=polyfun_with_QCed_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/%x%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename_qc.gz'
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/bellenguez'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final.tsv.gz'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/wightman_2021/wightman_2021_fixed.parquet'
#output='/gpfs/commons/home/tlin/output/wightman/wightman_all'

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'
brain_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K4me3/merged_annotations_ukb/brain_H3K4me3_seq_chr'
brain_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K27ac/merged_annotations_ukb/brain_H3K27ac_seq_chr'

#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
#output='/gpfs/commons/home/tlin/output/kunkle_all_2/all_anno'


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo


## create parquet file

#python munge_polyfun_sumstats.py \
#  --sumstats $summary_stats \
#  --out ~/data/bellenguez_munged.parquet 

summary_stats='/gpfs/commons/home/tlin/data/bellenguez_munged.parquet'

##1-2
#python polyfun.py \
#  --compute-h2-L2 \
#  --output-prefix $output \
#  --sumstats $summary_stats \
#  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno,$brain_H3K4me3,$brain_H3K27ac \
#  --w-ld-chr $bl/weights.UKB. \
#  --allow-missing

#echo finish polyfun1_2

#1-3
#for i in {1..22}
#do
#sbatch --export=chr=$i,output=$output /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh  
#done


#1-4
python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
   	--sumstats $summary_stats \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing

echo finish polyfun1_4
