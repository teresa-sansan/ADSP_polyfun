#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers

#file='Bellenguez_et_al_2021_hg37_no_dup.tsv.gz'
#qc_file_name='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_qc'

#file='Kunkle_et_al_2019_hg37_ldsc.tsv.gz'
#qc_file_name='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc'


#file='Wightman_et_al_2021_hg37_ldsc.tsv.gz'
file='processed/Wightman_2021_hg37_withbeta.tsv'
qc_file_name='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_withbeta'
ori=$(cat $file | wc -l)
echo There are $(expr $ori - 1) lines of SNPs in $file

## remove SNPs MAF < 0.01   ##kunkle and wightman don't have MAF info, so skip this step for kunkle
if false; then
echo 'filtering out MAF < 0.01...'
zcat $file | awk 'NR==1 || ($6>=0.01) {print}' > $qc_file_name.tsv
remove_maf=$(cat $qc_file_name.tsv| wc -l)
echo $(expr $ori - $remove_maf) 'of SNPs that MAF < 0.01 were removed, ' $(expr $remove_maf - 1) 'of SNPs remain.'
echo
fi


## duplicated SNPs   ##remember to change input

if true;then
echo "removing duplicated SNPs..."
cat $file | awk '{seen[$1]++; if(seen[$1]==1){print}}' > ${qc_file_name}_nodup.tsv
#cat $qc_file_name.tsv | awk '{seen[$1]++; if(seen[$1]==1){print}}' > ${qc_file_name}_nodup.tsv
nodup=$(cat ${qc_file_name}_nodup.tsv| wc -l)
echo $(expr $ori - $nodup) duplicated SNPs were removed,  $(expr $nodup - 1) of SNPs remain.
#echo $(expr $remove_maf - $nodup) duplicated SNPs were removed,  $(expr $nodup - 1) of SNPs remain.  ##if ran the first step, use this line

echo 
fi

## Ambiguous SNPs

echo "Removing ambiguous SNPs..."
cat ${qc_file_name}_nodup.tsv |awk '!(($4=="A" && $5 == "T") || ($4=="T" && $5 == "A") || ($4=="C" && $5 == "G") ||($4=="G" && $5 == "C")) {print}' > ${qc_file_name}.tsv 
no_ambiguous=$(cat ${qc_file_name}.tsv | wc -l)
echo $(expr $nodup - $no_ambiguous) 'ambiguous SNPs were removed'
echo total number of SNPs = $(expr $no_ambiguous - 1)


gzip $qc_file_name.tsv 
rm $qc_file_name_nodup.tsv
