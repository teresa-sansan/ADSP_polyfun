#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021
#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers
cd ~/data/qc_sumstat

file='Bellenguez_et_al_2021_hg37_no_dup.tsv.gz'
qc_file_name='Bellenguez_et_al_2021_hg37_qc'

ori=$(zcat $file | wc -l)
echo There are $(expr $ori - 1) lines of SNPs in $file

## remove SNPs MAF < 0.01
echo 'filtering out MAF < 0.01...'
zcat $file | awk 'NR==1 || ($6>=0.01) {print}' > $qc_file_name.tsv
remove_maf=$(cat $qc_file_name.tsv| wc -l)
echo $(expr $ori - $remove_maf) 'of SNPs that MAF < 0.01 were removed, ' $(expr $remove_maf - 1) 'of SNPs remain.'
echo


## duplicated SNPs
echo "removing duplicated SNPs..."
cat $qc_file_name.tsv | awk '{seen[$3]++; if(seen[$3]==1){print}}' > $qc_file_name_nodup.tsv
nodup=$(cat $qc_file_name_nodup.tsv| wc -l)
echo $(expr $remove_maf - $nodup) duplicated SNPs were removed,  $(expr $nodup - 1) of SNPs remain.
echo 

## Ambiguous SNPs

echo "Removing ambiguous SNPs..."
cat $qc_file_name_nodup.tsv |awk '!(($4=="A" && $5 == "T") || ($4=="T" && $5 == "A") || ($4=="C" && $5 == "G") ||($4=="G" && $5 == "C")) {print}' > $qc_file_name.tsv 
no_ambiguous=$(cat $qc_file_name.tsv | wc -l)
echo $(expr $nodup - $no_ambiguous) 'ambiguous SNPs were removed'
echo total number of SNPs = $(expr $no_ambiguous - 1)


gzip $qc_file_name.tsv 
rm $qc_file_name_nodup.tsv
