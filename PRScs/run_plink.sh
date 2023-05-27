#!/bin/bash
#SBATCH --job-name=prs_36k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/plink_output/log%x%j.log

#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/'
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/'
#cat $path/_pst_eff_a1_b0.5_phi1e-02_chr* > $path/agg_prscs_beta.txt



  for chr in {1..2}
  do

    chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    echo $chunk_num
    for ((i=1; i<=$chunk_num; i++)); 
    do
        ~/plink2/plink2 \
        --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/ADSP.chr${chr}.chunk${i} \
        --rm-dup 'force-first'
        --out $path/prscs/plink_output/chr${chr}.chunk${i}.dup

    done
  done






#  for chr in {1..5}
#  do
#     chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
#     echo $chunk_num
#     for ((i=1; i<=$chunk_num; i++)); 
#     do
#         ~/plink \
#         --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/ADSP.chr${chr}.chunk${i} \
#         --score $path/prscs/beta/chr${chr}.beta 2 4 6 \
#         --exclude $path/duplist.txt \
#         --q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
#         --out $path/prscs/plink_output/chr${chr}.chunk${i}

#         # for thres in e-5 0.001 0.005 0.01 0.05 0.1 0.5
#         # do
#         #     cat $path/prscs/plink_output/chr${chr}.chunk${i}.${thres}.profile  |tr -s ' '| cut -d ' ' -f 2-7 > $path/prscs/plink_output/chr${chr}.chunk${i}_${thres}.prs.tsv
#         # done
#     done
# done

## compare with directly running with plink (without PRSCS)
if false; then
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
--score  $path/agg_extract_0.3.tsv  2 4 12 \
--q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp \
--out /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/prs_plink

fi

# /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc
# /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc
# /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/agg_prscs_beta.txt 
# --score /gpfs/commons/home/tlin/output/wightman/prscs/all_anno/agg_prscs_beta.txt 2 4 6 \