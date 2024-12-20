#!/bin/bash
#SBATCH --job-name=finemap_no_partition_susie
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0824/susie/finemap/%x_%j.log

## double check if im running susie
cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}_*.npz"

## choose the sum stat you are running
##kunkle
if false; then
sumstat_name='kunkle'
sumstat="/gpfs/commons/home/tlin/output/kunkle/new_anno/${anno}/${anno}.${chr}.snpvar_ridge_constrained.gz"
n=63926
output="/gpfs/commons/home/tlin/output/kunkle/new_anno/${anno}/finemap"
fi

##bellenguez
if false; then
sumstat_name='bellenguez'
anno_path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
n=487511

##susie
#output="/gpfs/commons/home/tlin/output/bellenguez/new_sep22/susie/finemap"
fi

## wightman
if true; then
sumstat_name='wightman'
n=762971
anno_path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'

fi

## jansen
if false; then
sumstat_name='jansen'
sumstat='/gpfs/commons/home/tlin/output/jansen/jansen'
n=450734
#output='/gpfs/commons/home/tlin/output/jansen/finemap'
output='/gpfs/commons/home/tlin/output/jansen/susie'
fi

echo run ${sumstat_name}

for i in $FILES
do	
	ld=$(echo $i | cut -d'.' -f 1 ) 
	filename=$(echo $i | cut -d'/' -f 10 |cut -d '.' -f 1)
	start=$(echo $filename| cut -d '_' -f 2)
	end=$(echo $filename| cut -d '_' -f 3)
 ## snpvar_constrained.gz
	if true; then
		if [ $max_num_snp -eq 1 ] 
		then
		python finemapper.py \
					--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_ridge.gz \
					--n $n --chr ${chr} --start $start --end $end \
					--method susie \
					--max-num-causal 1 \
					--allow-missing \
					--non-funct \
					--out $anno_path/susie/finemap/max_snp_${max_num_snp}/${sumstat_name}.${chr}.$start.$end.gz
		else
		python finemapper.py \
			--ld $ld \
			--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_ridge.gz \
			--n $n --chr ${chr} --start $start --end $end \
			--method susie \
			--max-num-causal $max_num_snp \
			--allow-missing \
			--non-funct \
			--out $anno_path/susie/finemap/max_snp_${max_num_snp}/${sumstat_name}.${chr}.$start.$end.gz
		fi
	fi

done