cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed
sumstat='bellenguez_hg38_parquet_flipped.tsv'

for chr in {1..22}
do
echo running chr$chr ...
echo -e "ID\tCHR\tPOS\tRef\tAlt\tSNP\tN\tMAF\tZ\tPval" > carma/chr${chr}.tsv
tail -n+2 $sumstat | awk -v chr=$chr '($1 == chr){print $1 ":" $2 ":" $4 ":" $5, $1, $2, $4, $5, $3, $9, $10, $12, $8}' OFS="\t" >> carma/chr${chr}.tsv
gzip carma/chr${chr}.tsv

done