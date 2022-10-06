#!/bin/bash
#SBATCH --job-name=bellenguez_rescue
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=50:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_sep22/bl/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld"
sumstat="/gpfs/commons/home/tlin/output/bellenguez/new_sep22/$anno/"
missing="/gpfs/commons/home/tlin/output/bellenguez/new_sep22/bl/finemap/max_snp_${max_snp}/missing_window.txt"

#sed -i s/'chr'/''/g /gpfs/commons/home/tlin/output/bellenguez/new_sep22/bl/finemap/max_snp_10/missing_window.txt
#sed -i s/'chr'/''/g /gpfs/commons/home/tlin/output/bellenguez/new_sep22/bl/finemap/max_snp_5/missing_window.txt
## polyfun, not susie
## remember to change output name

for line in $(cat $missing)
do	 
	if  grep -q ".gz" <<< $line ; then
		chr=$(echo $line| awk -F '.'  '{print $2}')
		start=$(echo $line| awk -F '.' '{print $3}')
		end=$(echo $line|awk -F '.' '{print $4}')
	    ldfile=$(echo $line| awk -F '.' '{print "chr" $2 "_" $3 "_" $4}')

   		echo finemapping on chr${chr} BP ${start} to ${end}
		python finemapper.py \
			--sumstats $sumstat/${anno}.${chr}.snpvar_constrained.gz \
			--n 487511 \
			--ld $FILES/${ldfile} \
	  		--chr ${chr} --start $start --end $end \
	  		--method susie \
     	  	--max-num-causal ${max_snp} \
	  		--allow-missing \
			--out $sumstat/finemap/max_snp_${max_snp}/${anno}.chr${chr}.${start}.${end}.gz
	fi
	

done


