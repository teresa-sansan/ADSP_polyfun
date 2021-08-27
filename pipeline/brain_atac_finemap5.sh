#!/bin/bash
#SBATCH --job-name=brain_atac_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bl_brain_atac/finemap/correct_max_causal_5/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"

for i in $FILES
do	
	
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	ld_ref=$(echo $i|cut -d'.' -f1)

	python finemapper.py \
		--ld $ld_ref \
		--sumstats /gpfs/commons/home/tlin/output/bl_brain_atac/kunkle_brain_atac.${chr}.snpvar_ridge_constrained.gz \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 5 \
	  	--allow-missing \
		--out /gpfs/commons/home/tlin/output/bl_brain_atac/finemap/correct_max_causal_5/finemap_bl_brainatac.${chr}.$start.$end.gz \
	

done
