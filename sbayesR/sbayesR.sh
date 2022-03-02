#!/bin/bash
#SBATCH --job-name=sbayesR
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/sbayesR/kunkle/fixed_0224/%x_%j.log
#SBATCH --cpus-per-task=20


for chr in 1
do
/gpfs/commons/home/tlin/sbayesR/gctb_2.03beta_Linux/gctb  --sbayes R \
     --ldm /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/SBayesR/ld_matrix/band_ukb_10k_hm3/band_chr${chr}.ldm.sparse \
     --gwas-summary /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_no_dup.ma \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --maf 0.001 \
     --exclude-mhc \
     --thread 20 \
     --hsq 0.5 \
     --chain-length 35000 \
     --burn-in 1000 \
     --seed 9448 \
     --out /gpfs/commons/home/tlin/output/sbayesR/test_chr${chr}
done

