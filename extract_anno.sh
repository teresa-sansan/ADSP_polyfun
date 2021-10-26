#!/bin/bash
#SBATCH --job-name=extract_anno
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/%x%j.log 

cd ~/polyfun_omer_repo/
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
path='/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1'
prefix="all_anno"

chr=22
##annotations
baseline='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB'
all='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/combined_AD_annotations_polyfun/combined_AD_annotations_polyfun_' 
brain_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K4me3/merged_annotations_ukb/brain_H3K4me3_seq_chr'
brain_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_H3K27ac/merged_annotations_ukb/brain_H3K27ac_seq_chr'


echo $baseline.${chr}.annot.parquet
for i in $FILES
do 
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
        start=$(echo $filename| cut -d'_' -f 2)
        end=$(echo $filename| cut -d'_' -f 3)

	python extract_annotations.py \
	--pips $path/${prefix}.${chr}.${start}.${end}.gz \
	--annot ${baseline}.${chr}.annot.parquet \
	--pip-cutoff 0.5 \
	--allow-missing \
	--out $path/${prefix}_extract_${chr}.$start.$end.csv
	
	break

done   




#	--annot /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.${chr}.annot.parquet,/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr${chr}.annot.gz \
# ${all}${chr}.annot.gz, ${brain_H3K4me3}${chr}.annot.gz, ${brain_H3K27ac}${chr}.annot.gz \
