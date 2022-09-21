##0920_22
## Chirag create new bellenguez file
## But that file has "chr" in the chr column, which polyfun can't recognize
## So this script is just (1) create a copy to my repo (2) remove 'chr' from chr (e.g. change chr1 to 1)
 
cp /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_et_al_2021_hg37.tsv.gz /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_et_al_2021_hg37_new_sep20.tsv.gz
mv /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_et_al_2021_hg37_new_sep20.tsv.gz 
zcat Bellenguez_et_al_2021_hg37_new_sep20.tsv.gz| awk  -v FS='\t' -v OFS='\t' '{gsub (/'chr'/,"",$2); print$0'} > Bellenguez_et_al_2021_hg37_new_sep20_fix_chr.tsv /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/
