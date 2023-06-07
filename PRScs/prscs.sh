#!/bin/bash
#SBATCH --job-name=wightman_prscs_17k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer/%x_%j.log


cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun



# python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
#     --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
#     --sst_file=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_4prscs_not0.tsv \
#     --n_gwas=762971 --phi=1e-2 \
#     --out_dir=/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/PIP_not0/ 


python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
    --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
    --sst_file=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/finemap/max_snp_5/agg_4prscs.tsv \
    --n_gwas=762971 --phi=1e-2 \
    --out_dir=/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer/ 



# path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
# chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
# echo chumk num is $chunk_num


# for ((i=1; i<=$chunk_num; i++)); 
# do
# echo $path/ADSP.chr${chr}.chunk${i} 

# ~/plink \
# python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
#    --bim_prefix=$path/plink_hg38/ADSP.chr${chr}.chunk${i} \
#    --sst_file=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv \
#    --chrom=$chr --n_gwas=762971 --phi=1e-2 --out_dir=$path/PRS_hg38/prscs/chr${chr}_chunk${i}
# done

# cat $path/PRS/prscs/chr${chr}_chunk* > $path/PRS/prscs/chr${chr}.beta 
# echo "concat all chunks in chr " $chr






## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv
## --ref_dir: LD panel file
## bim_prefix: bim file for target data
## sst_file: sumstat_file 

## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc