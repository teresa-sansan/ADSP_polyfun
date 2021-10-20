#!/bin/bash
#SBATCH --job-name=bellenguez_rescue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun


FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld"
sumstat='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl'
missing=$sumstat/finemap_susie/max_snp_$max_snp/missing_window.txt


##susie!

for line in $(cat $missing)
do	 
	if  grep -q ".gz" <<< $line ; then
		chr=$(echo $line| awk -F '.'  '{print $2}')
		start=$(echo $line| awk -F '.' '{print $3}')
		end=$(echo $line|awk -F '.' '{print $4}')
	        ldfile=$(echo $line| awk -F '.' '{print "chr" $2 "_" $3 "_" $4}')

   		echo finemapping on chr${chr} BP ${start} to ${end}
		python finemapper.py \
			--sumstats $sumstat/bellenguez_bl.${chr}.snpvar_constrained.gz \
			--n 487511 \
			#--ld $FILES/${ldfile} \
			--non-func \
	  		--chr ${chr} --start $start --end $end \
	  		--method susie \
     	  		--max-num-causal 1 \
	  		--allow-missing \
			--out $sumstat/finemap_susie/max_snp_1/finemap_bellenguez_susie.${chr}.$start.$end.gz
	fi
	

done

#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all"
#missing="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${max_snp}/missing_window.txt"



