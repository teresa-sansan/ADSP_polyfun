touch IBSS_not_converge.txt
#finemap_path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'
finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_update/finemap'

cd $finemap_path
pwd
touch IBSS_not_converge.txt

for max_snp in 1 3 5 7 10
do
echo max_snp $max_snp | tee -a $finemap_path/IBSS_not_converge.txt 	
	for i in {1..22} 
	do 
		converge=$(cat max_snp_$max_snp/bellenguez.chr$i.*.gz.log| grep converge| wc -l)
		if [ $converge != 0 ]
		then
			echo chr $i, not converged = $converge | tee -a $finemap_path/IBSS_not_converge.txt
		fi
  
	done 

done

