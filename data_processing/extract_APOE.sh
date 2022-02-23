cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021
ori_sum_stat='bellenguez_2021_final_rename.tsv'
qc_sum_stat='bellenguez_2021_final_qc.gz'


#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019
#ori_sum_stat='kunkle_etal_Stage1_info.tsv'
#qc_sum_stat='Kunkle_etal_Stage1_qc.gz'



cat $ori_sum_stat |head -1 > APOE_only_no_qc.tsv
cat $ori_sum_stat |awk '{if($1==19 && $2>45409011 && $2<45412650) print $0}'>> APOE_only_no_qc.tsv

zcat $qc_sum_stat |head -1 > APOE_only_qc.tsv
zcat $qc_sum_stat |awk '{if($1==19 && $2>45409011 && $2<45412650) print $0}'>> APOE_only_qc.tsv


##or
# zcat bellenguez_2021_final.tsv.gz |grep chr19 |awk '{if($2>45409011 && $2<45412650) print $0}' 
