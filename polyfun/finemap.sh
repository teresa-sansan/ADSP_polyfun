#!/bin/bash
#SBATCH --job-name=finemap_bellenguez_new_anno
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=40G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/only_ml/finemap/%x_%j.log

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
if true; then
sumstat_name='bellenguez'
anno_path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/new_anno/${anno}/${anno}.${chr}.snpvar_ridge_constrained.gz"
#output="/gpfs/commons/home/tlin/output/bellenguez/new_anno/${anno}/finemap"
n=487511

##susie
#output="/gpfs/commons/home/tlin/output/bellenguez/new_sep22/susie/finemap"
fi

## wightman
if false; then
sumstat_name='wightman'
n=762971
#sumstat="/gpfs/commons/home/tlin/output/wightman/wightman_all.${chr}.snpvar_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/bl/bl.${chr}.snpvar_constrained.gz"
anno_path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'
#sumstat="/gpfs/commons/home/tlin/output/wightman/new_anno_0203/${anno}/${anno}.${chr}.snpvar_ridge_constrained.gz"
#output='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/all/finemap'
#output="/gpfs/commons/home/tlin/output/wightman/new_anno_0203/${anno}/finemap"
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
	
	
	if false; then
		for anno in bl bl_dl_annotations bl_brain_atac
		do	
			## --sumstats ${sumstat}.${chr}.snpvar_constrained.gz
			## $output/max_snp_${max_num_snp}/${sumstat_name}.chr${chr}.$start.$end.gz
		if [ $max_num_snp -eq 1 ] 
		then
		python finemapper.py \
					--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
					--n $n \
					--chr ${chr} --start $start --end $end \
					--method susie \
					--max-num-causal 1 \
					--allow-missing \
					--out $output/${anno}/max_snp_${max_num_snp}/{anno}.chr${chr}.$start.$end.gz
		else
		python finemapper.py \
			--ld $ld \
			--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
			--n $n \
			--chr ${chr} --start $start --end $end \
			--method susie \
			--max-num-causal $max_num_snp \
			--allow-missing \
			--out $output/${anno}/max_snp_${max_num_snp}/${anno}.chr${chr}.$start.$end.gz
		fi
		done
	fi

	if true; then
		if [ $max_num_snp -eq 1 ] 
		then
		python finemapper.py \
					--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
					--n $n --chr ${chr} --start $start --end $end \
					--method susie \
					--max-num-causal 1 \
					--allow-missing \
					--out $anno_path/$anno/finemap/max_snp_${max_num_snp}/${sumstat_name}.${chr}.$start.$end.gz
		else
		python finemapper.py \
			--ld $ld \
			--sumstats $anno_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
			--n $n --chr ${chr} --start $start --end $end \
			--method susie \
			--max-num-causal $max_num_snp \
			--allow-missing \
			--out $anno_path/$anno/finemap/max_snp_${max_num_snp}/${sumstat_name}.${chr}.$start.$end.gz
		fi
	fi

done