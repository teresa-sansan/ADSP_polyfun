#!/bin/bash
#SBATCH --job-name=extract_anno
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --output=output/bl_roadmap_microglia/extract_anno/test/%x%j.log 

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
PIPFILES="bl_roadmap_microglia"

for i in $FILES
do 
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
        start=$(echo $filename| cut -d'_' -f 2)
        end=$(echo $filename| cut -d'_' -f 3)

	python extract_annotations.py \
	--pips /gpfs/commons/home/tlin/polyfun/output/$PIPFILES/finemap/finemap_$PIPFILES.${chr}.$start.$end.gz \
	--annot /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.${chr}.annot.parquet \
	--pip-cutoff 0.5 \
	--allow-missing \
	--out /gpfs/commons/home/tlin/polyfun/output/$PIPFILES/extract_anno/test/bl_anno_${chr}_$start.$end.csv

done   




#	--annot /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.${chr}.annot.parquet,/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr${chr}.annot.gz \
