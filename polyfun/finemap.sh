#!/bin/bash
#SBATCH --job-name=finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=30:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/finemap_fixed_assertion_susie_iter/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*.gz"

##kunkle
if false; then
echo run kunkle
sumstat="/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/kunkle.${chr}.snpvar_constrained.gz"
n=63926
output='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_${max_num_snp}/kunkle'
fi


##bellenguez
if true; then
echo run bellenguez
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.${chr}.snpvar_constrained.gz"
n=487511
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap_fixed_assertion_susie_iter/'
fi

## wightman
if false; then
echo run wightman
sumstat="/gpfs/commons/home/tlin/output/wightman/wightman_all.${chr}.snpvar_constrained.gz"
n=74004
output="/gpfs/commons/home/tlin/output/wightman/finemap_fixed_assertion_susie_iter/"
fi

for i in $FILES
do	
	filename=$(echo $i | cut -d'/' -f 10 |cut -d '.' -f 1)
	start=$(echo $filename| cut -d '_' -f 2)
	end=$(echo $filename| cut -d '_' -f 3)
	
	if [ $max_num_snp -eq 1 ] 
	then
	python finemapper.py \
                --sumstats $sumstat \
                --n $n \
                --chr ${chr} --start $start --end $end \
                --method susie \
                --max-num-causal $max_num_snp \
                --allow-missing \
                --out $output/max_snp_${max_num_snp}/chr${chr}.$start.$end.gz
	else
	python finemapper.py \
		--ld $filename \
		--sumstats $sumstat \
		--n $n \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal $max_num_snp \
	  	--allow-missing \
		--out $output/max_snp_${max_num_snp}/chr${chr}.$start.$end.gz 
	fi
done
