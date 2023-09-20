#!/bin/bash
#SBATCH --job-name=prscs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=25G
#SBATCH --time=3:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSCS/36k_ibd_adsp_fixed/%x_%j.log

'''
recommend to run it per chr chunk. (each chr takes ~30 min). run it with run_prscs.sh directly.
'''

cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

sumstat_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed'


python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
    --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr} \
    --sst_file=$sumstat_path/$sumstat \
    --chr $chr \
    --n_gwas=$n_sumstat --phi=1e-2 \
    --out_dir=$output_dir/


## old geno
## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc


## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv
## --ref_dir: LD panel file
## bim_prefix: bim file for target data
## sst_file: sumstat_file 