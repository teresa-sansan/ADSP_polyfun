#!/bin/bash
#SBATCH --job-name=finemap_fixed_convergence_issue_rerun
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=100:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld"

##bellenguez
if true; then
echo run bellenguez
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/bellenguez"
n=487511
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/'
fi


## run it using 1 MB window sliding window, with 0.5 MB overlap.
## this file only have the failing regions in max_num_snp  = 10
## set max_num_snp to 3 (10/3 = 3)
max_num_snp=3

while IFS= read -r line
do	
	
	chr=$(echo $line | cut -d '.' -f 1|sed s/finemap_max_snp_3_chr//)
	start=$(echo $line | cut -d '.' -f 2)
	end=$(echo $line | cut -d '.' -f 3)
	
	test_500000=$(expr $start % 1000000)

	
	if [[ test_500000 -eq 1 ]] ; then
		LD_start=$(expr $start - 1000000)
		LD_end=$(expr $start + 2000000)
	else
		echo not full
		LD_start=$(expr $start - 1500000)
		LD_end=$(expr $start + 1500000)
	fi
	

	#block_tail=$(echo $line | cut -d '.' -f 2| cut -d '_' -f 2)   ##only taking the tail 1MB, here is only taking the end of the 1MB tail
	#block_head=$(expr $block_tail - 1000000) 
	#block_1=$(expr $block_head - 500000)
	#block_2=$(expr $block_head + 500000)
	#block_3=$(expr $block_tail + 500000) 
	#LD_start=$(expr $block_head - 1000000)
	#LD_end=$(expr $block_head + 2000000)
		
	if [ $end > 0 ]; then
		python finemapper_max_iter_1000.py \
		--ld $FILES/chr${chr}_${LD_start}_${LD_end} \
		--sumstats $sumstat.${chr}.snpvar_constrained.gz \
		--n $n 	--chr $chr --start $start --end $end \
		--method susie \
    		--max-num-causal ${max_num_snp} \
		--allow-missing \
		--out $output/try_rescue_not_converge/finemap_max_snp_${max_num_snp}_chr${chr}.${start}.${end}.gz 
	fi
	
done < /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/no_result.txt
