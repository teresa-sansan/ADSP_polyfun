#!/bin/bash
#SBATCH --job-name=bellenguez_convergence
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=100:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_sep22/susie/finemap/max_snp_10/try_rescue_not_converge/%x_%j.log


## change #56 and add the revision in #57.
cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld"

anno='all_anno'
##bellenguez
if true; then
echo run bellenguez
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/bellenguez"
sumstat='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/'
n=487511
output='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/susie/finemap/max_snp_10/rescue_not_convergence/'
echo "run not converge regions in" ${path} ${anno}'/max_snp_10/IBSS_not_converge_list.txt'

#			--sumstats $sumstat/$anno/${anno}.${chr}.snpvar_constrained.gz \
#			--out $output/$anno/max_snp_10/try_rescue_not_converge/finemap_max_snp_${max_num_snp}_chr${chr}.${start}.${end}.gz 

fi


#wightman
if false; then
echo run wightman
sumstat='/gpfs/commons/home/tlin/output/wightman/fixed_0224'
n=74004
output='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10'
echo "run not converge regions in /gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/run_IBSS_not_converge_list.txt"

#			--sumstats $sumstat/wightman_all.${chr}.snpvar_constrained.gz \

fi


## run it using 1 MB window sliding window, with 0.5 MB overlap.
## this file only have the failing regions in max_num_snp  = 10
## set max_num_snp to 3 (10/3 = 3)
max_num_snp=3

while IFS= read -r line
do	
	
	chr=$(echo $line | cut -d '.' -f 1|sed s/chr//)
	#block_tail=$(echo $line | cut -d '.' -f 2| cut -d '_' -f 2)   ##only taking the tail 1MB, here is only taking the end of the 1MB tail
	block_tail=$(echo $line | cut -d '.' -f 3| cut -d ' ' -f 1)   ##only taking the tail 1MB, here is only taking the end of the 1MB tail, re-modified

	block_head=$(expr $block_tail - 1000000) 
	block_1=$(expr $block_head - 500000)
	block_2=$(expr $block_head + 500000)
	block_3=$(expr $block_tail + 500000) 
	LD_start=$(expr $block_head - 1000000)
	LD_end=$(expr $block_head + 2000000)
		
	for start in  $block_1 $block_head $block_2   ## do finemap in each chunk. 
	do
	end=$(expr $start + 1000000)
		if [ $start > 0 -a ! -f $output/$anno/finemap/max_snp_10/try_rescue_not_converge/bellenguez_max_snp_${max_num_snp}.chr${chr}.${start}.${end}.log ]; then
			python finemapper_max_iter_1000.py \
			--ld $FILES/chr${chr}_${LD_start}_${LD_end} \
			--sumstats $sumstat/$anno/${anno}.${chr}.snpvar_constrained.gz \
			--n $n 	--chr $chr --start $start --end $end \
	  		--method susie \
			--non-funct \
    	  	--max-num-causal ${max_num_snp} \
	  		--allow-missing \
			--out $sumstat/susie/finemap/max_snp_10/try_rescue_not_converge/finemap_max_snp_${max_num_snp}.chr${chr}.${start}.${end}.gz 

		fi
	
	done
done < $sumstat/susie/finemap/max_snp_10/run_IBSS_not_converge_list.txt

#/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/run_IBSS_not_converge_list.txt
