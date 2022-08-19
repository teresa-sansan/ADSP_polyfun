## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

for chr in 19
do
~/plink \
--bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_${chr} \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_APOE_qc.tsv 1 4 6 \
--q-score-range range_list.txt $sumstats_path/aggregate.pvalue  \
--extract /gpfs/commons/home/tlin/output/cT/new_plink/wightman/fixed_0224/qc/chr${chr}.valid.snp \
--out /gpfs/commons/home/tlin/output/cT/genomewide_plink/kunkle/APOE_qc/kunkle_APOE_qc_chr${chr}
done

#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/APOE_2SNPs_qc.tsv 3 4 6 header \
#--out /gpfs/commons/home/tlin/output/kunkle/kunkle_qc/2SNP_qc.prs
