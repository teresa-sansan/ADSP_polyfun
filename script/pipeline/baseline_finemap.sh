#!/bin/bash
#SBATCH --job-name=baseline_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=output/kunkle_UKBBbaseline/output/finemap/%x%j.log

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
	#echo $start 
	#echo $end
	#echo ${chr}
	
	python finemapper.py \
		--sumstats output/kunkle_UKBBbaseline/output/kunkle_UKBBbaseline.${chr}.snpvar_constrained.gz \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 5 \
	  	--allow-missing \
		--out ~/polyfun/output/kunkle_UKBBbaseline/output/finemap/max_causal_5/finemap_UKBBbaseline.${chr}.$start.$end.gz \
	

done
