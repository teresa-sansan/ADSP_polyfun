#!/bin/bash
#SBATCH --job-name=susie_wightman
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=15:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/susie/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}_*.npz"

##kunkle
if false; then
echo run kunkle
sumstat="/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl/bl.${chr}.snpvar_constrained.gz"
n=63926
output='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl'
fi


##bellenguez
if false; then
echo run bellenguez
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl/bl"
n=487511
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/'
fi

## wightman
if true; then
echo run wightman
sumstat="/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/bl/bl"
n=74004
output="/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/susie"
fi

for i in $FILES
do	
	ld=$(echo $i | cut -d'.' -f 1 ) 
	filename=$(echo $i | cut -d'/' -f 10 |cut -d '.' -f 1)
	start=$(echo $filename| cut -d '_' -f 2)
	end=$(echo $filename| cut -d '_' -f 3)

	#for anno in bl bl_dl_annotations bl_brain_atac
	#do
	anno='susie'
	if [ $max_num_snp -eq 1 ] 
	then
	python finemapper.py \
                --sumstats ${sumstat}.${chr}.snpvar_constrained.gz \
                --n $n \
                --chr ${chr} --start $start --end $end \
                --method susie \
                --non-funct \
                --max-num-causal $max_num_snp \
                --allow-missing \
                --out  $output/polypred/max_snp_${max_num_snp}_susie_chr${chr}.$start.$end.gz
	else
	python finemapper.py \
		--ld $ld \
		--sumstats ${sumstat}.${chr}.snpvar_constrained.gz \
		--n $n \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
		--non-funct \
     	  	--max-num-causal $max_num_snp \
	  	--allow-missing \
		--out $output/polypred/max_snp_${max_num_snp}_susie_chr${chr}.$start.$end.gz 
	fi
	#done
done
