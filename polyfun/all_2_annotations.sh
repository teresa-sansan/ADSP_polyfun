#!/bin/bash
#SBATCH --job-name=new_bellenguez
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=180G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_sep22/%x%j.log

#sumstat_name='wightman'
#sumstat_name='jansen'
sumstat_name='bellenguez/new_sep22'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_et_al_2021_hg37.tsv.gz'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Jansen_et_al_2019_hg37_ldsc.tsv.gz'

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'
brain_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K4me3/merged_annotations_ukb/brain_H3K4me3_seq_chr'
brain_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K27ac/merged_annotations_ukb/brain_H3K27ac_seq_chr'


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo

## create parquet file

if false; then
#  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_previous.munged.parquet
python munge_polyfun_sumstats.py \
  --sumstats $summary_stats \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstat_name}_et_al_2021_hg37_ldsc.munged.parquet
fi
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman_previous.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/kunkle_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Wightman.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/jansen_et_al_2021_hg37_ldsc.munged.parquet'
summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sep22_new_bellenguez_et_al_2021_hg37_ldsc.munged.parquet'


#output='/gpfs/commons/home/tlin/output/jansen'
output='/gpfs/commons/home/tlin/output/bellenguez'


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
echo finish polyfun1_2
fi


##1-2_anno
if false; then
output='bl/bl'
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix /gpfs/commons/home/tlin/output/$sumstat_name/$output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB. \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

output='bl_dl_annotations/bl_dl_annotations'
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix /gpfs/commons/home/tlin/output/$sumstat_name/$output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$all_anno \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

output='bl_brain_atac/bl_brain_atac'
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix /gpfs/commons/home/tlin/output/$sumstat_name/$output \
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
for i in {11..22}
do                                                                                                                                                                                                                                                                         
  #sbatch --export=output=/gpfs/commons/home/tlin/output/$sumstat_name/fixed_0224_annotations/$anno/$anno /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh
  sbatch --export=chr=$i,output=/gpfs/commons/home/tlin/output/$sumstat_name/$anno/$anno /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh
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
for i in bl bl_dl_annotations bl_brain_atac bl
do
python polyfun.py \
        --compute-h2-bins \
        --output-prefix /gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/$i/$i \
        --sumstats $summary_stats \
        --w-ld-chr $bl/weights.UKB. \
        --allow-missing

echo finish polyfun1_4
done
fi