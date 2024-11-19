#!/bin/bash
#SBATCH --job-name=prs_36k_adsp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=3:00:00
#SBATCH --array=1-22%18
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/prscs/prs_middlefile//%x%j.log

note=$(cat << EOF
running adsp need to use hg38!
running ukbb need to use hg19!
EOF
)

echo "$note"
chr=$SLURM_ARRAY_TASK_ID
#chr=$1
## ADSP LD
#path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_adsp_ld_panel/bellenguez'
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/prscs'
#bellenguez_snp_hg38 ='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/bellenguez_hg38_qc.pvalue'
sumstat_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/'
bellenguez_snp_hg38='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/prscs/bellenguez_hg38_parquet_flipped_4prscs.snp'

~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.chr${chr} \
--score $path/bellenguez_pst_eff_a1_b0.5_phi1e-02_chr${chr}.txt 2 4 6 \
--out $path/prs_middlefile/chr${chr}.qc
#--q-score-range ../plink/range_list.txt $sumstat_path/bellenguez_hg38_parquet_flipped_a1a2.snp \
#--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt $bellenguez_snp_hg38 \


## UKBB LD
# echo "use UKBB LD"
# path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez_rerun_0909/'
# bellenguez_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc.pvalue'
# #wightman_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp'
# ~/plink \
# --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr} \
# --score $path/bellenguez_pst_eff_a1_b0.5_phi1e-02_chr${chr}.txt 2 4 6 \
# --q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt $bellenguez_snp \
# --out $path/prs_middlefile/chr${chr}.qc


for thres in e-6 e-5 e-4 0.001 0.01 0.1 0.9
#for chr in {1..22}
do
    cat $path/prs_middlefile/chr${chr}.qc.${thres}.profile  |tr -s ' '| cut -d ' ' -f 2-7 > $path/prs_middlefile/chr${chr}.qc_${thres}.prs.tsv      
    #cat $path/prs_middlefile/chr${chr}.qc.profile  |tr -s ' '| cut -d ' ' -f 2-7 > $path/prs_middlefile/chr${chr}.qc.prs.tsv      
    #rm $path/prs_middlefile/chr${chr}.qc.${thres}.profile 
done


# ## compare with directly running with plink (without PRSCS)
# if false; then
# ~/plink \
# --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
# --score  $path/agg_extract_0.3.tsv  2 4 12 \
# --q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
# --out /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prs_plink

# fi
# --exclude /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/duplist.txt \
# /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc
# /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr${chr}

# --score /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/agg_prscs_beta.txt 2 4 6 \

#--exclude $path/duplist.txt --extract $path/vcf_rsid/uniq_ADSP.chr${chr}.chunk${i}.snp \

##ukb_ld:/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr}