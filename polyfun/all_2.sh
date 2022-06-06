#!/bin/bash
#SBATCH --job-name=polyfun_bellenguez_update0224
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/%x%j.log


summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_ldsc.tsv.gz'
##create munge
#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_test/bellenguez'
#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez'
#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/bellenguez'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_no_dup.tsv.gz'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename_qc.gz'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final.tsv.gz'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Wightman_et_al_2021_hg37_ldsc.tsv.gz'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/wightman_2021/wightman_2021_fixed.parquet'
#output='/gpfs/commons/home/tlin/output/wightman/previous/wightman_all'

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

if false;then
#  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_previous.munged.parquet


python munge_polyfun_sumstats.py \
  --sumstats $summary_stats \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_ldsc.munged.parquet

fi
summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_previous.munged.parquet'

##1-2
if false; then
echo $summary_stats

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno,$brain_H3K4me3,$brain_H3K27ac \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

#echo finish polyfun1_2

fi



##1-2_anno
if false; then
#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl'
#python polyfun.py \
#  --compute-h2-L2 \
#  --output-prefix $output \
#  --sumstats $summary_stats \
#  --ref-ld-chr $bl/baselineLF2.2.UKB. \
#  --w-ld-chr $bl/weights.UKB. \
#  --allow-missing

#output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl_dl_annotations'
#python polyfun.py \
#  --compute-h2-L2 \
#  --output-prefix $output \
#  --sumstats $summary_stats \
#  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno \
#  --w-ld-chr $bl/weights.UKB. \
#  --allow-missing

output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl_brain_atac/bl_brain_atac'
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$brain_H3K4me3,$brain_H3K27ac \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

fi



#1-3
if false; then
for i in {1..22}
do
sbatch --export=chr=$i,output=$output /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh  
done
fi


#1-3 anntations
if false; then   
for anno in bl bl_dl_annotations bl_brain_atac
do                                                                                                                                                                                                                                                                          
  for i in {1..22}
  do
  sbatch --export=chr=$i,output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/$anno/$anno /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh
  done
done
fi

#1-4
if false; then
python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
   	--sumstats $summary_stats \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing

echo finish polyfun1_4
fi


if true; then
for i in bl bl_dl_annotations bl_brain_atac
do
python polyfun.py \
        --compute-h2-bins \
        --output-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/$i/$i \
        --sumstats $summary_stats \
        --w-ld-chr $bl/weights.UKB. \
        --allow-missing

echo finish polyfun1_4
done
fi

