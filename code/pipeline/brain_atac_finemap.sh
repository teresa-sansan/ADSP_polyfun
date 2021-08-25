#!/bin/bash
#SBATCH --job-name=brain_atac_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/polyfun/output/kunkle_UKBB_brain_atac/finemap/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"

for i in $FILES
do	
	#echo $i
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	python finemapper.py \
		--sumstats output/kunkle_UKBB_brain_atac/kunkle_brain_atac.${chr}.snpvar_ridge_constrained.gz \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 1 \
	  	--allow-missing \
		--out /gpfs/commons/home/tlin/polyfun/output/kunkle_UKBB_brain_atac/finemap/finemap_UKBB_brainatac.${chr}.$start.$end.gz \
	

done
