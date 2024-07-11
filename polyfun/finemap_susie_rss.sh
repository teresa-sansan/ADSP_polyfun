#!/bin/bash
#SBATCH --job-name=finemap_susie_RSS_test
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/%x_%j.log


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

chr=19
echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}_*.npz"

## choose the sum stat you are running

##bellenguez hg38
sumstat_name='bellenguez_hg38'
sumstat='bellenguez_hg38_qc.tsv'

##bellenguez hg19

anno_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/bellenguez'

n=487511
prefix='bellenguez'
anno='baseline'
bp=44893972


# start=$((bp - 100000))
# end=$((bp + 100000))

echo run functional informed finemapping
python /gpfs/commons/home/tlin/script/polyfun/finemapper_susie_rss.py \
					--geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
					--sumstats ${anno_path}_${anno}/${prefix}_${anno}.${chr}.snpvar_ridge.gz \
					--n $n --chr ${chr} --start $((bp - 100000)) --end $((bp + 100000)) \
					--method susie \
					--max-num-causal 10 \
					--allow-missing \
					--debug-dir /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/ \
					--out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/test.txt


