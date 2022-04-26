## remove the two most significant APOE allele: rs42938 and rs7412

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed

cat Kunkle_et_al_2019_hg37_ldsc.tsv | awk '{if ($2 =='19' && $3 >= 45409011 && $3 <= 45412650) print $0}' > Kunkle_APOE.tsv 
cat Kunkle_et_al_2019_hg37_ldsc_qc.tsv | awk '{if ($2 =='19' && $3 >= 45409011 && $3 <= 45412650) print $0}' > Kunkle_APOE_qc.tsv 


grep -Fvxf Kunkle_APOE.tsv Kunkle_et_al_2019_hg37_ldsc.tsv > Kunkle_remove_APOE.tsv
grep -Fvxf Kunkle_APOE_qc.tsv Kunkle_et_al_2019_hg37_ldsc_qc.tsv > Kunkle_remove_APOE_qc.tsv

