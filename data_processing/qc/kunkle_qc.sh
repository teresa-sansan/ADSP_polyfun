cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019

## duplicated SNPs
cat Kunkle_etal_Stage1_results.txt | awk '{seen[$3]++; if(seen[$3]==1){print}}'|gzip > Kunkle_etal_Stage1_results_nodup.gz


## Ambiguous SNPs
zcat Kunkle_etal_Stage1_results_nodup.gz |awk '!(($4=="A" && $5 == "T") || \
 ($4=="T" && $5 == "A") || ($4=="C" && $5 == "G") ||($4=="G" && $5 == "C")) {print}' |gzip > Kunkle_etal_Stage1_qc.gz 

echo CHR	BP	SNP	A1	A2	Beta	SE	Pvalue > Kunkle_etal_Stage1_qc_renamed.tsv
zcat Kunkle_etal_Stage1_qc.gz | tail -n+2|tr ' ' '\t' >> Kunkle_etal_Stage1_qc_renamed.tsv && gzip Kunkle_etal_Stage1_qc_renamed.tsv
rm  Kunkle_etal_Stage1_qc.gz
mv Kunkle_etal_Stage1_qc_renamed.tsv.gz  Kunkle_etal_Stage1_qc.gz
