#!/bin/bash
#SBATCH --job-name=Finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/%x%j.log

cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr1*" 
summary_stat="/gpfs/commons/home/tlin/output/kunkle_all_2/all_anno.1.snpvar_constrained.gz"

start=190000001

while [ $start -le 249000001 ]
do	
	
	end=`expr $start + 3000000`
	echo start=$start, end=$end

	python finemapper.py \
		--sumstats $summary_stat \
		--n 63926 \
	  	--chr 1 --start $start --end $end \
	  	--method susie \
		--max-num-causal 1 \
	  	--allow-missing \
		--out /gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/all_anno.$1.$start.$end.gz 
	start=`expr $start + 1000000`	

done
