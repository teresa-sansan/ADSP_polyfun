cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021

## duplicated SNPs
echo "removing duplicated SNPs..."
ori=$(zcat bellenguez_2021_final_rename.tsv.gz | wc -l)
zcat bellenguez_2021_final_rename.tsv.gz | awk '{seen[$3]++; if(seen[$3]==1){print}}' > bellenguez_2021_nodup.tsv
nodup=$(cat bellenguez_2021_nodup.tsv| wc -l)
echo $(expr $ori - $nodup) duplicated SNPs were removed
echo 

## Ambiguous SNPs

echo "Removing ambiguous SNPs..."
cat bellenguez_2021_nodup.tsv |awk '!(($4=="A" && $5 == "T") || \
 ($4=="T" && $5 == "A") || ($4=="C" && $5 == "G") ||($4=="G" && $5 == "C")) {print}' |gzip > bellenguez_2021_final_rename_qc.gz 

gzip bellenguez_2021_nodup.tsv
no_ambiguous=$(zcat bellenguez_2021_final_rename_qc.gz | wc -l )
echo $(expr $nodup - $no_ambiguous) ambiguous SNPs were removed
echo total number of SNPs = $(expr $no_ambiguous - 1)
