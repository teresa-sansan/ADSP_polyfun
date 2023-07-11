#!/bin/bash
#SBATCH --job-name=prs_36k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez/%x%j.log

#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/'
#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/'
path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez'
#touch  $path/agg_prscs_beta.txt

#cat $path/_pst_eff_a1_b0.5_phi1e-02_chr* > $path/agg_prscs_beta.txt
wightman_snp=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp
bellenguez_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc.pvalue'


 for chr in {21..22}
 
 do
        ~/plink \
        --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/ADSP_qc_chr${chr} \
        --score $path/chr${chr}_pst_eff_a1_b0.5_phi1e-02_chr${chr}.txt 2 4 6 \
        --exclude /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/duplist.txt \
        --q-score-range /gpfs/commons/home/tlin/script/plink/range_list.txt $bellenguez_snp  \
        --out $path/chr${chr}.qc

        for thres in e-5 0.001 0.005 0.01 0.05 0.1 0.5 
        do
            cat $path/chr${chr}.qc.${thres}.profile  |tr -s ' '| cut -d ' ' -f 2-7 > $path/chr${chr}.qc_${thres}.prs.tsv
            
        done
done

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

#--exclude $path/duplist.txt --extract $path/vcf_rsid/uniq_ADSP.chr${chr}.chunk${i}.snp \