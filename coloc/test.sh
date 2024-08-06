#!/bin/bash
#SBATCH --job-name=finemap_eqtl
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/%x_%j.log


source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun


eqtl_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene'
gene='ENSG00000159140.21'
chr=21
bp=34567305
start=$((bp - 1000000))
end=$((bp + 1000000))
python  ~/polyfun_omer_repo/finemapper.py \
        --sumstats $eqtl_path/${gene}_sumstats.hg38.parquet \
        --geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr${chr} \
        --n 346 --chr ${chr} --start $start --end $end \
        --method susie --max-num-causal 10 \
        --non-funct --allow-missing \
	--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/${gene}.rds \
        --susie-resvar 1 --out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/${gene}.txt
