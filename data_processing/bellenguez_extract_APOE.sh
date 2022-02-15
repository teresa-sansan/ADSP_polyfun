cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021

zcat bellenguez_2021_final_rename.tsv.gz |head -1 > bellenguez_2021_fin_APOE.tsv
zcat bellenguez_2021_final_rename.tsv.gz |awk '{if($1==19 && $2>45409011 && $2<45412650) print $0}'>> bellenguez_2021_fin_APOE.tsv


##or
# zcat bellenguez_2021_final.tsv.gz |grep chr19 |awk '{if($2>45409011 && $2<45412650) print $0}' 
