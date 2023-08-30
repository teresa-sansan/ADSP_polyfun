#!/bin/bash 
## run  bellenguez n fixed 

for i in {1..20}
do
	for max_snp in 1 5 
	do	    
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"
        	#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='bl' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='all_except_enformer' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='no_ml_new' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='all_anno' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='all_enformer' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='no_ml' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='bl' finemap.sh
			
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='all' finemap.sh
			#sbatch --export=chr=$i,max_num_snp=${max_snp},anno='no_ml' finemap.sh
			sbatch --export=chr=$i,max_num_snp=${max_snp},anno='only_ml' finemap.sh

	done
done
