#!/bin/bash
#SBATCH --job-name=wightman_prscs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=25:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/%x%j.log


cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun


python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
    --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
    --sst_file=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_4prscs.tsv \
    --n_gwas=762971 --phi=1e-2 --out_dir=/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/



## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv
## --ref_dir: LD panel file
## bim_prefix: bim file for target data
## sst_file: sumstat_file 

