cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019
tail Kunkle_etal_Stage1_results.txt -n+2 |cut -d ' ' -f 3,8| grep -v "NA" > Kunkle_no_qc.SNP  
