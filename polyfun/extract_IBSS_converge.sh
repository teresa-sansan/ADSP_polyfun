## Teresa Lin on 04_22_22
## purpose: SuSiE's IBSS algorithm if not converged, will generate result that is not trust-worthy. 
## 	    So runs a job to check how many blocks fail to converge is essential. 
## note: remember to check finemap file name before running this. 


#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl_dl_annotations'
#finemap_path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/old/finemap'
finemap_path='/gpfs/commons/home/tlin/output/jansen/susie'
#finemap_path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/all_anno/finemap'
#finemap_path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/'

# summary_stat='susie'
#summary_stat='bellenguez'
summary_stat='jansen'
cd $finemap_path
echo block counts for those failed to converge > $finemap_path/IBSS_not_converge_count.txt

for max_snp in 10
do
#echo WRITE not converge list in max_snp = $max_snp | tee -a $finemap_path/IBSS_not_converge_count.txt
echo max_snp $max_snp | tee  $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt 	
	for chr in {1..22} 
	do	
		dont_converge='F'
		for i in $finemap_path/max_snp_$max_snp/${summary_stat}.chr${chr}.*.gz.log
		do 
		converge=$(cat $i| grep converge| wc -l)
		if [ $converge != 0 ]
		then
			ld=$( echo $i| cut -d '/' -f 10)   ## sometmes it's 12 bellenguez ## 10 for jansen ## 13 for wightman
			pos=$( echo $ld| cut -d '.' -f 2-4)
			echo $pos ';' $converge 'jobs fail to converge.'| tee -a $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt
			echo ' ' | tee -a $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt
			dont_converge='T'
		fi
		done
		if [ $dont_converge == 'T' ]
		then
		count=$(cat $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt | grep chr$chr| wc -l)
		echo chr$chr : $count >> $finemap_path/IBSS_not_converge_count.txt

		fi
		
	done
	echo ' ' >> IBSS_not_converge_count.txt
done

cat $finemap_path/max_snp_$max_snp/IBSS_not_converge_list.txt| grep chr > $finemap_path/max_snp_$max_snp/run_IBSS_not_converge_list.txt


