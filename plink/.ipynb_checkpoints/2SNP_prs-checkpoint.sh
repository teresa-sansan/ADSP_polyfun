~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/19_filt \
--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove2SNP_noqc.tsv 1 4 6 header \
--out /gpfs/commons/home/tlin/output/kunkle/kunkle_qc/updated_0224_remove_2SNP_noqc.prs


#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_remove2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_2SNP_qc.tsv 1 4 6 header \
#--score /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/APOE_2SNPs_qc.tsv 3 4 6 header \
#--out /gpfs/commons/home/tlin/output/kunkle/kunkle_qc/2SNP_qc.prs
