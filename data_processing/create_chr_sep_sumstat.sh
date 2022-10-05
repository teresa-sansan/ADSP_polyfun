
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed 
file='wightman_fixed_beta_qc.tsv'
save_file='wightman_chr_sep/wightman_fixed_beta_qc_chr'
for chr in {1..4}
do
echo start chr ${chr}
head -1 $file > ${save_file}${chr}.tsv
cat $file | awk -v chr="$chr" ' $2 == chr' >> ${save_file}${chr}.tsv
cat ${save_file}${chr}.tsv|cut -f 1,7 > ${save_file}${chr}.pvalue
done

