#!/bin/bash
#SBATCH --job-name=assertionerror_tryout
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=42:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/%x%j.log
cd ~/polyfun_omer_repo

#--job-name=susie_kunkle

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
anno="bl_brain_atac"

#chr=10
#max_num_snp=10
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*.npz" 
summary_stat="/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl_brain_atac/bl_brain_atac.${chr}.snpvar_constrained.gz"
output="/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations"

for i in $FILES
do	
	
	ld=$(echo $i | cut -d'.' -f 1 )
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)

 
	if [ $max_num_snp -eq 1 ]
	then
		python re-clone/finemapper.py \
		--sumstats $summary_stat \
		--n 63926 \
	  	--chr $chr --start $start --end $end \
	  	--method susie \
		--max-num-causal $max_num_snp \
		--non-funct \
	  	--allow-missing \
		--out $output/$anno/max_snp_${max_num_snp}/chr${chr}.$start.$end.gz 
	else	
		python re-clone/finemapper.py \
                --ld $ld\
                --sumstats $summary_stat \
                --n 63926 \
                --chr $chr --start $start --end $end \
                --method susie \
                --max-num-causal $max_num_snp \
		--non-funct \
                --allow-missing \
		--out $output/$anno/max_snp_${max_num_snp}/chr${chr}.$start.$end.gz 
	fi	

done
