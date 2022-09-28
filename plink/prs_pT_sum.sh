
#file_name='ADSP_no_apoe_qc_all'
#dir_path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/kunkle/ADSP_no_apoe'

dir_path='/gpfs/commons/home/tlin/output/cT/new_plink_genomewide/jansen/ADSP_qc_all/'
file_name='ADSP_qc_all'

#for dir in bellenguez kunkle wightman
#do
	#cd /gpfs/commons/home/tlin/output/cT/new_plink/$dir/fixed_0224/polyfun_beta_no_clump
	#cd $path/$qc
 	#cd  /gpfs/commons/home/tlin/output/cT/genomewide_plink/$dir/$qc
	cd $dir_path
	pwd
	for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
	do
	echo write pT_$i.prs
	awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' ${file_name}_chr*.$i.profile > pT_${i}.prs
	done
	echo 
 #done
#done
echo finish summing up all PRS!
