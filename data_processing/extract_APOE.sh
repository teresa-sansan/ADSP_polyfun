#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021
#ori_sum_stat='bellenguez_2021_final_rename.tsv'
#qc_sum_stat='bellenguez_2021_final_qc.gz'


#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019
#ori_sum_stat='kunkle_etal_Stage1_info.tsv'
#qc_sum_stat='Kunkle_etal_Stage1_qc.gz'


cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers
sum_stat='Bellenguez_et_al_2021_hg37.tsv.gz'

zcat $sum_stat |head -1 > Bellenguez_APOE_only_no_qc.tsv
zcat Bellenguez_et_al_2021_hg37.tsv.gz |awk '{if($2=="chr19" && $3>45409011 && $3<45412650) print $0}'>>Bellenguez_APOE_only_no_qc.tsv

echo ttest

if(false);then

cat $ori_sum_stat |head -1 > APOE_only_no_qc.tsv
cat $ori_sum_stat |awk '{if($1==19 && $2>45409011 && $2<45412650) print $0}>>APOE_only_no_qc.tsv

zcat $qc_sum_stat |head -1 > APOE_only_qc.tsv
zcat $qc_sum_stat |awk '{if($1==19 && $2>45409011 && $2<45412650) print $0}'>> APOE_only_qc.tsv
fi

##or
# zcat bellenguez_2021_final.tsv.gz |grep chr19 |awk '{if($2>45409011 && $2<45412650) print $0}' 
