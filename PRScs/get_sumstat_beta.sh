wightman='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta_qc.tsv'
snp='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_4prscs_not0.tsv'

awk 'NR==FNR{a[$1]=$0; next} $1 in a {print $1,$4,$5,$10,$7}' $snp $wightman > /gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/prscs_beta_from_sumstat.tsv