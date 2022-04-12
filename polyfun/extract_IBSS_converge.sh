touch IBSS_not_converge.txt
#finemap_path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'
finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_updated/finemap'
#summary_stat='wightman'
summary_stat='bellenguez'
cd $finemap_path
echo check convergence | tee $finemap_path/IBSS_not_converge_list.txt

for max_snp in 1 3 5 7 10
do
echo max_snp $max_snp | tee -a $finemap_path/IBSS_not_converge_list.txt 	
	for chr in {1..22} 
	do
		for i in $finemap_path/max_snp_$max_snp/${summary_stat}.chr$chr.*.gz.log
		do 
		converge=$(cat $i| grep converge| wc -l)
		if [ $converge != 0 ]
		then
			#echo $converge 'fails to converge.'| tee -a  $finemap_path/IBSS_not_converge_list.txt
			#echo ' ' | tee -a $finemap_path/IBSS_not_converge_list.txt
			ld=$( echo $i| cut -d '/' -f 11)
			pos=$( echo $ld| cut -d '.' -f 2-3)
			echo $pos | tee -a $finemap_path/IBSS_not_converge_list.txt
		fi
		done
	done
echo ' ' | tee -a $finemap_path/IBSS_not_converge_list.txt 

done
