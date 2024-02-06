#!/bin/bash
#SBATCH --job-name=wightman_all
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G 
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0824_no_partitions/bl/%x%j.log

## --partition bigmem

#output='/gpfs/commons/home/tlin/output/kunkle/new_anno/all_anno/all_anno'
#output='/gpfs/commons/home/tlin/output/kunkle/new_anno/all_enformer/all_enformer'
#output='/gpfs/commons/home/tlin/output/kunkle/new_anno/no_ml/no_ml'
#output='/gpfs/commons/home/tlin/output/kunkle/new_anno/bl/bl'

#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/update_all+enformer' 
#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/all_except_enformer' 
#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/bl/bl' 
#output='/gpfs/commons/home/tlin/output/bellenguez/new_anno/bl/bl'

output='/gpfs/commons/home/tlin/output/wightman/new_anno_0824_no_partitions/'

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
# all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'

## new_anno_0203
# deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
# enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/enformer/enformer_high_h2_chr'
# glass_lab='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab/glass_lab_high_h2_chr'
# glass_lab_enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer/glass_lab_enformer_high_h2_chr'
# roadmap='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/roadmap/roadmap_high_h2_chr'

## new_anno_0822
# bl_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/baseline/baseline_high_h2_chr'
# deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
# glass_lab_enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer/glass_lab_enformer_high_h2_chr'
# roadmap='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/roadmap/roadmap_high_h2_chr'

## new_anno0824
bl_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/baseline/baseline_high_h2_chr'
deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/enformer/enformer_high_h2_chr'
glasslab='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab/glass_lab_high_h2_chr'
glass_lab_enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer/glass_lab_enformer_high_h2_chr'
roadmap='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/roadmap/roadmap_high_h2_chr'

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo

## create parquet file

if false; then
python munge_polyfun_sumstats.py \
  --sumstats $summary_stats \
  --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstat_name}_fixed_beta.munged.parquet
fi
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/kunkle_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/jansen_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sep22_new_bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.munged.parquet'

##1-2
if true; then
echo running $summary_stats

python polyfun.py \
  --compute-h2-L2 \
  --no-partitions \
  --output-prefix $output/no_ml/no_ml \
  --sumstats $summary_stats \
  --ref-ld-chr $bl_anno,$roadmap,$glasslab \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing &

python polyfun.py \
  --compute-h2-L2 \
  --no-partitions \
  --output-prefix $output/only_ml/only_ml \
  --sumstats $summary_stats \
  --ref-ld-chr $bl_anno,$deepsea,$enformer,$glass_lab_enformer \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing &

python polyfun.py \
  --compute-h2-L2 \
  --no-partitions \
  --output-prefix $output/all/all \
  --sumstats $summary_stats \
  --ref-ld-chr $bl_anno,$roadmap,$glasslab,$deepsea,$enformer,$glass_lab_enformer \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing 



echo finish polyfun1_2
fi

#1-3
if false; then
for i in {1..22}   
do
sbatch --export=chr=$i,output=$output/bl/bl /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh  
done
fi

#1-4
if false; then
python polyfun_assertion_error.py \
  --compute-h2-bins \
  --output-prefix $output/bl/bl \
  --sumstats $summary_stats \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

echo finish polyfun1_4
fi



 