## chr9,12,18,22 have missing file in p<5e-8

#cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT
#cd /gpfs/commons/home/tlin/output/cT/kunkle/qc_check
#cd /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc

#path='/gpfs/commons/home/tlin/output/cT/wightman/'
#path='/gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224'
#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224'
#path='/gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224/qc_check_target'

#for qc in before_qc qc_on_base qc_on_target qc
#for qc in qc_all_maf01 qc_target_maf01	

#for qc in qc_on_individual qc_on_variant
for path in /gpfs/commons/home/tlin/output/cT/wightman /gpfs/commons/home/tlin/output/cT/bellenguez/fixed_0224 /gpfs/commons/home/tlin/output/cT/kunkle/fixed_0224
do	
	#echo in $qc
	#cd $path/$qc
	cd $path/qc_on_variant_sumstat
	pwd
	for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
	do
	echo write pT_$i.prs
	awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' *.$i.profile > pT_${i}.prs
	done
	echo '  '
done

echo finish summing up all PRS!
