#!/bin/bash
#SBATCH --job-name=susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=stdout_%j.log


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

cd ~/polyfun

python finemapper.py \
  --sumstats data/kunkle_2019_teresa/AD_Kunkle_etal_Stage1.parquet \
  --n 63926\
  --chr 1 \
  --start $start \
  --end $end \
  --method susie \
  --max-num-causal 1 \
  --non-funct \
  --allow-missing \
  --out ~/polyfun/output/UKBBoutput/ukbb_chr1_$start.$end.gz
