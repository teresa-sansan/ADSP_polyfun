#!/bin/bash
#SBATCH --job-name=polyloc_stage2_
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G 
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/polyloc/%x%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

polyfun_path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/'

## Stage 1
# python polyloc_reclone/polyloc.py \
#     --compute-partitions \
#     --output-prefix $polyfun_path/polyloc/polyloc \
#     --posterior $polyfun_path/max_snp_5/aggregate.all.txt \
#     --bfile-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/UKB_snp_list/UKB_SNP_list_chr
    
## Stage 2
python polyloc_reclone/polyloc.py \
    --output-prefix $polyfun_path/polyloc/polyloc \
    --compute-ldscores \
    --ld-ukb \
    --ld-dir /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/ \
    --chr $chr


   # 
## add this so it wont take that much time. If dont add this tag, the default one is 29
# 36k /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr 

# 19k /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_chr
# --bfile-chr /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr \
# --num-bins 20 
#  --bfile-chr /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/UKB_snp_list/UKB_SNP_list_chr 