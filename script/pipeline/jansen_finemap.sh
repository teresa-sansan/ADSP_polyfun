#!/bin/bash
#SBATCH --job-name=jansen_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --output=output/jansenetal/%x%j.log

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
	#echo $start 
	#echo $end
	#echo ${chr}
	
	python finemapper.py \
		--sumstats output/jansenetal/bl_jansen.${chr}.snpvar_constrained.gz   \
		--n 450734 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 1 \
	  	--out ~/polyfun/output/jansenetal/finemap/finemap_jansen.${chr}.$start.$end.gz \
		--allow-missing
	

done
