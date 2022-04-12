#!/bin/bash
#SBATCH --job-name=finemap_fixed_convergence_issue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=30:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/finemap/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/"

##bellenguez
if true; then
echo run bellenguez
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/bellenguez"
n=487511
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/finemap/'
fi


max_num_snp=10

for line in $(cat /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/finemap/IBSS_not_converge_list.txt | tail -n 168)
do	
	chr=$(echo $line | cut -d '.' -f 1|sed s/chr//)
	LD=$(echo $line | cut -d '.' -f 2| cut -d '_' -f 2)
	start=$(expr $LD - 1000000) 
	end=$(expr $LD + 2000000)
	echo $chr, $LD
	echo chr${chr}_${start}_$end
	python finemapper.py \
		--ld $FILES/chr${chr}_${start}_$end \
		--sumstats $sumstat.${chr}.snpvar_constrained.gz \
		--n $n \
	  	--chr $chr --start $start --end $LD \
	  	--method susie \
     	  	--max-num-causal $max_num_snp \
	  	--allow-missing \
		--out $output/max_snp_${max_num_snp}/test_convergence_chr${chr}.$start.$LD.gz 
done
