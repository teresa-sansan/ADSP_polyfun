#!/bin/bash
#SBATCH --job-name=Polypred_mergebeta
#SBATCH --mail-type=FAIL,end
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=60G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/w_prscs/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo/
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

output_path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polypred/w_prscs/'

## step1, calculate merge beta
# python polypred.py --combine-betas \
#   --betas /gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/aggregate.all.txt,/gpfs/commons/home/tlin/output/wightman/prscs/original/agg_beta_header.txt \
#   --pheno /gpfs/commons/home/tlin/data/17K_pheno.tsv \
#   --output-prefix $output_path/merged \
#   --plink-exe ~/plink \
#   /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc.bed

#  --betas /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz,~/output/sbayesR/bellenguez_agg_snpRes \
#  --betas ~/output/sbayesR/bellenguez_agg_snpRes,/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/aggregrate.all.txt.gz \


## step2, calculate PRS
python polypred.py \
    --predict \
    --betas $output_path/merged.betas \
    --output-prefix $output_path/polypred.predictions \
    --plink-exe ~/plink \
    /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc.bed
