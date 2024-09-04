#!/bin/bash
#SBATCH --job-name=prscs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=25:00:00
#SBATCH --array=1-20%20
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSCS/36k_adsp_ld_panel/bellenguez/%x_%j.log


##recommend to run it per chr chunk. (each chr takes ~30 min). run it with run_prscs.sh directly.


cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

chr=$SLURM_ARRAY_TASK_ID

n_bellenguez=487511
sumstat_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed'
sumstat='prscs/bellenguez_hg38_parquet_flipped_4prscs.tsv'
#sumstat='prscs/bellenguez_hg38_chr21.tsv'

n_sumstat=487511
output_dir='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_adsp_ld_panel/bellenguez'

python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_ADSP36K_4PRScs_Aug/ldblk_adsp_chr \
    --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.chr${chr} \
    --sst_file=$sumstat_path/$sumstat \
    --chr $chr --n_gwas=$n_sumstat --phi=1e-2 \
    --out_dir=$output_dir


## ukb_ld /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr}
## ukb_ref /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur

## old geno
## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc


## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv
## --ref_dir: LD panel file
## bim_prefix: bim file for target data
## sst_file: sumstat_file 