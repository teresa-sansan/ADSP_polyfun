#!/bin/bash
#SBATCH --job-name=Finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle_all/finemap/%x%j.log

cd ~/polyfun

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*" 
summary_stat="/gpfs/commons/home/tlin/output/kunkle_all/all_anno_.${chr}.snpvar_constrained.gz"


for i in $FILES
do	

	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	python finemapper.py \
		--sumstats $summary_stat \
		--n 63926 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
		--max-num-causal 1 \
	  	--allow-missing \
		--out /gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno.${chr}.$start.$end.gz 
	

done
