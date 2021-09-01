#!/bin/bash
#SBATCH --job-name=baseline_finemap_5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=output/bl/output/finemap/max_causal_5/%x_%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
for i in $FILES
do	
	#echo $i
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	echo $filename
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	ld_ref=$(echo $i|cut -d'.' -f1)
	
	python finemapper.py \
		--ld $ld_ref \
		--sumstats output/bl/output/kunkle_UKBBbaseline.${chr}.snpvar_constrained.gz \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 5 \
	  	--allow-missing \
		--out ~/polyfun/output/bl/output/finemap/max_causal_5/finemap_bl.${chr}.$start.$end.gz \
	

done
