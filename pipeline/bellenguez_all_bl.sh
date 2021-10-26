#!/bin/bash
#SBATCH --job-name=bellenguez
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/%x%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/bellenguez_bl'

#summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
#output='/gpfs/commons/home/tlin/output/kunkle_all_2/all_anno'


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
cd ~/polyfun_omer_repo

##1-2
python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB. \
  --w-ld-chr $bl/weights.UKB. \
  --allow-missing

#echo finish polyfun1_2

##1-3
for i in {1..22}
do
  python polyfun.py \
      --compute-ldscores \
      --output-prefix $output \
      --ld-ukb \
      --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld
      --chr $i
done
echo finish polyfun1_3

##1-4
python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
    	--sumstats $summary_stats \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing

echo finish polyfun1_4
