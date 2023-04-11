#!/bin/bash 
## run  wightman n fixed 

for i in {6..22}
do
	for max_snp in  1 5 10
	do	
		echo "run_max_snp_$max_snp"
        	echo "run chr$i"

		if true; then
        	sbatch --export=chr=$i,max_num_snp=${max_snp},anno='glasslab' finemap.sh
        fi
	done
done
