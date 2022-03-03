## remove the two most significant APOE allele: rs42938 and rs7412

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/
zcat Kunkle_et_al_2019_hg37_ldsc.tsv.gz | awk '{if ($1 != "rs429358" && $1 != "rs7412") print $0}' > processed/Kunkle_remove2SNP_noqc.tsv



cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed
zcat Kunkle_et_al_2019_hg37_ldsc_qc.tsv.gz | awk '{if ($1 != "rs429358" && $1 != "rs7412") print $0}' > Kunkle_remove2SNP_qc.tsv 

