#!/bin/bash
#SBATCH --job-name=roadmap_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/specific_col/finemap_susie/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
sumstat="/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/specific_col/bl_roadmap.${chr}.snpvar_ridge_constrained.gz"

for i in $FILES
do	
	#echo $i
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	python finemapper.py \
		--sumstats $sumstat \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 1 \
	  	--allow-missing \
		--non-funct \
		--out "/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/specific_col/finemap_susie/finemap_bl_roadmap.${chr}.$start.$end.gz"

	

done
