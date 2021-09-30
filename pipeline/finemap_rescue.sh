#!/bin/bash
#SBATCH --job-name=bellenguez_rescue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun


FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld"
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all"
missing="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${max_snp}/missing_window.txt"

for line in $(cat $missing)
do	 
	if  grep -q ".gz" <<< $line ; then
		chr=$(echo $line| awk -F '.'  '{print $2}')
		start=$(echo $line| awk -F '.' '{print $3}')
		end=$(echo $line|awk -F '.' '{print $4}')
	        ldfile=$(echo $line| awk -F '.' '{print "chr" $2 "_" $3 "_" $4}')

   		echo finemapping on chr${chr} BP ${start} to ${end}
		python finemapper.py \
			--ld $FILES/${ldfile} \
			--sumstats $sumstat.${chr}.snpvar_constrained.gz \
			--n 487511 \
	  		--chr ${chr} --start $start --end $end \
	  		--method susie \
     	  		--max-num-causal ${max_snp} \
	  		--allow-missing \
			--out "/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${max_snp}/finemap_bellenguez_all_2.${chr}.$start.$end.gz"
	fi
	

done
