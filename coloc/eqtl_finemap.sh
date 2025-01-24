#!/bin/bash
#SBATCH --job-name=eqtl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/%x_%j.log

gene=$1
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
### make sure run this when susieR version is 0.11.92

#remove.packages("susieR")
#remotes::install_github("stephenslab/susieR@v0.11.92")

#finemap on eQTL
eqtl_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene'

# echo $output
# if [  -f $eqtl_path/${gene}_sumstats_hg38.gz ]; then
# 	echo "renaming header for eQTL ..."
# 	#awk 'BEGIN{OFS=FS="\t"} {sub(/^chr/, "", $3)} 1' $eqtl_path/microglia_eqtl_chr${chr}.tsv > file_chr${chr}.tmp
#     zcat $eqtl_path/${gene}_sumstats_hg38.gz | sed -e '1s/chr/CHR/' -e '1s/pos/BP/'  -e '1s/variant_id/SNP/' -e '1s/ref/A1/'  -e '1s/alt/A2/'  -e '1s/fixed_z/Z/' > $eqtl_path/${gene}_sumstats_hg38
# else
#     echo ${gene}_sumstats_hg38.gz not found
# fi

gene_count='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene_count.tsv'
n=$(awk -F '\t' -v gene="$gene" '$1 == gene {print $2}' "$gene_count")
info=$(python get_parquet_tss.py $gene)
read chr tss <<< "$info"
echo $chr
window=1000000

output='chr'${chr}'_'${gene}
echo $output

start_try=$((tss - window))
if [ "$start_try" -lt 0 ]; then
    start=0
else
    start=$start_try
fi

end=$((tss + window))

python /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_susie_rss.py \
					--sumstats $eqtl_path/${gene}_sumstats.hg38.parquet \
                    --geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
					--n $n --chr ${chr} --start $start --end $end \
					--method susie \
					--max-num-causal 10 \
					--non-funct \
					--susie-resvar 1 \
					--no-sort-pip \
					--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_$output \
					--allow-missing \
					--out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_${output}.txt

