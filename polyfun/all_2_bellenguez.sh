#!/bin/bash
#SBATCH --job-name=bellenguez_enformer
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=300G 
#SBATCH --time=15:00:00
#SBATCH --partition bigmem
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno/all_enformer/%x%j.log

sumstat_name='bellenguez'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_ldsc.tsv.gz'
##create munge
#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/update_all+enformer' 
#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/all_except_enformer' 

#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/bl/bl' 
output='/gpfs/commons/home/tlin/output/bellenguez/new_anno/all_enformer/all_enformer'

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
all_anno='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_'

## new_anno_0203
deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/enformer/enformer_high_h2_chr'
glass_lab='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab/glass_lab_high_h2_chr'
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
summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sep22_new_bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.munged.parquet'

##1-2
if false; then
echo running $summary_stats

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$enformer,$glass_lab_enformer \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing
echo finish polyfun1_2
fi

#1-3
if false; then
for i in {1..22}   
do
sbatch --export=chr=$i,output=$output /gpfs/commons/home/tlin/script/polyfun/polyfun_1_3.sh  
done
fi

#1-4
if true; then
python polyfun_assertion_error.py \
 --compute-h2-bins \
  --output-prefix $output \
  --sumstats $summary_stats \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

echo finish polyfun1_4
fi


