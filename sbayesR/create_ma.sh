cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers
echo "SNP	A1	A2	freq	b	se	p	N" > Bellenguez_et_al_2021_hg37_no_dup.ma
zcat Bellenguez_et_al_2021_hg37_no_dup.tsv.gz| tail -n +2 | cut -f 1,4,5,6,7,8,9,14 >>Bellenguez_et_al_2021_hg37_no_dup.ma
