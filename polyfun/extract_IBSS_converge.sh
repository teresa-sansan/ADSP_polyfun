## Teresa Lin on 04_22_22
## purpose: SuSiE's IBSS algorithm if not converged, will generate result that is not trust-worthy. 
## 	    So runs a job to check how many blocks fail to converge is essential. 
## note: remember to check finemap file name before running this. 



#finemap_path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'
finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/old/finemap'
#summary_stat='wightman'
summary_stat='bellenguez'

cd $finemap_path
echo block counts for those failed to converge > IBSS_not_converge_count.txt

for max_snp in  3 5 7 10
do
echo WRITE not converge list in max_snp = $max_snp | tee -a IBSS_not_converge_count.txt
echo max_snp $max_snp | tee  $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt 	
	for chr in {1..22} 
	do	
		dont_converge='F'
		for i in $finemap_path/max_snp_$max_snp/finemap_${summary_stat}.$chr.*.gz.log
		do 
		converge=$(cat $i| grep converge| wc -l)
		if [ $converge != 0 ]
		then
			ld=$( echo $i| cut -d '/' -f 12)
			pos=$( echo $ld| cut -d '.' -f 2-4)
			echo chr$pos ';' $converge 'jobs fail to converge.'| tee -a $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt
			echo ' ' | tee -a $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt
			dont_converge='T'
		fi
		done
		if [ $dont_converge == 'T' ]
		then
		count=$(cat $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt | grep chr$chr| wc -l)
		echo chr$chr : $count >> IBSS_not_converge_count.txt

		fi
		
	done
	echo ' ' >> IBSS_not_converge_count.txt
done



