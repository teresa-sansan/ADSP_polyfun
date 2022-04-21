cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021
#ori_sum_stat='bellenguez_2021_final_rename.tsv'
qc_sum_stat='bellenguez_2021_final_qc.tsv'

##most centric
#awk 'FNR==NR {lines[$3]; next} $2 in lines {print $0}' $qc_sum_stat /gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP.tsv > Bellenguez_qc_interested_SNP.tsv

##MAX
awk 'FNR==NR {lines[$3]; next} $2 in lines {print $0}' $qc_sum_stat /gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP_maxPIP.tsv > Bellenguez_qc_interested_SNP_maxpip.tsv

##min
awk 'FNR==NR {lines[$3]; next} $2 in lines {print $0}' $qc_sum_stat /gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP_minPIP.tsv > Bellenguez_qc_interested_SNP_minPIP.tsv

#awk 'FNR==NR {lines[$3]; next} $2 in lines {print $0}' $ori_sum_stat /gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP.tsv > Bellenguez_interested_SNP.tsv
