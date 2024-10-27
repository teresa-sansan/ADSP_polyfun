#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/%x_%j.log


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
### make sure run this when susieR version is 0.11.92

#remove.packages("susieR")
#remotes::install_github("stephenslab/susieR@v0.11.92")

chr=$1
bp=$2
anno=$3

output='chr'${chr}'_'${bp}'_'$anno
#echo $output
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}_*.npz"

## choose the sum stat you are running

##bellenguez hg38
sumstat_name='bellenguez_hg38'
sumstat='bellenguez_hg38_qc.tsv'

#anno_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_backup_teresa/bellenguez'
anno_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/teresa/backup_oct/'

n=487511
prefix='bellenguez'


start=$((bp - 1000000))
end=$((bp + 1000000))


if [ $anno == 'susie' ]; then
	echo run finemapping on susie
	python /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_susie_rss.py \
						--geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
						--sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/bellenguez_hg38_parquet_flipped_a1a2.tsv \
						--n $n --chr ${chr} --start $start --end $end \
						--method susie \
						--max-num-causal 10 \
						--non-funct \
						--no-sort-pip \
						--susie-resvar 1 \
						--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/${anno}_$output \
						--allow-missing \
						--out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/${anno}_${output}.tsv


else
	echo run functional informed finemapping, on $anno
	python /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_susie_rss.py \
						--geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
						--sumstats ${anno_path}/bellenguez_${anno}.${chr}.snpvar_ridge.gz \
						--n $n --chr ${chr} --start $start --end $end \
						--method susie \
						--max-num-causal 10 \
						--no-sort-pip \
						--susie-resvar 1 \
						--allow-missing \
						--out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/${anno}_${output}.tsv \
						--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/snp_of_interest/${anno}_${output} \

fi



#finemap on eQTL
# eqtl_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene'
# gene=$3

# output='chr'${chr}'_'${bp}'_'${gene}
# echo $output
# # if [  -f $eqtl_path/${gene}_sumstats_hg38.gz ]; then
# # 	echo "renaming header for eQTL ..."
# # 	#awk 'BEGIN{OFS=FS="\t"} {sub(/^chr/, "", $3)} 1' $eqtl_path/microglia_eqtl_chr${chr}.tsv > file_chr${chr}.tmp
# #     zcat $eqtl_path/${gene}_sumstats_hg38.gz | sed -e '1s/chr/CHR/' -e '1s/pos/BP/'  -e '1s/variant_id/SNP/' -e '1s/ref/A1/'  -e '1s/alt/A2/'  -e '1s/fixed_z/Z/' > $eqtl_path/${gene}_sumstats_hg38
# # else
# #     echo ${gene}_sumstats_hg38.gz not found
# # fi


# python /gpfs/commons/home/tlin/polyfun_omer_repo//munge_eqtl_sumstats.revised.py \
#     --sumstats $eqtl_path/${gene}_sumstats_hg38.gz --n 346 --out $eqtl_path/${gene}_sumstats.hg38.parquet


# python  ~/polyfun_omer_repo/finemapper.py \
#         --sumstats $eqtl_path/${gene}_sumstats.hg38.parquet \
#         --geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
#         --n 346 --chr ${chr} --start $start --end $end \
#         --method susie --max-num-causal 10 \
#         --non-funct --allow-missing \
# 	--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/${gene}.rds \
#         --susie-resvar 1 --out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/${gene}.txt


# python /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_susie_rss.py \
# 					--sumstats $eqtl_path/${gene}_sumstats.hg38.parquet \
#                     --geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
# 					--n 346 --chr ${chr} --start $start --end $end \
# 					--method susie \
# 					--max-num-causal 10 \
# 					--non-funct \
# 					--susie-resvar 1 \
# 					--no-sort-pip \
# 					--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/eqtl_z2_$output \
# 					--allow-missing \
# 					--out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/eqtl_z2_${output}.txt
	

#--debug /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/${output}_ \
##--ld /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/susie_chr${chr}_${bp}_susie.ld \