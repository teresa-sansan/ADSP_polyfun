#path="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_susie"
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
prefix=bellenguez.chr
#prefix=finemap_bellenguez_all2_susie
#path="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie"
#prefix=finemap_bellenguez_susie
if [ -f $path/check_finemap.txt ]; then rm $path/check_finemap.txt; fi #remove pre-exist file
touch $path/check_finemap.txt   ##create file

function check
{
for i in {1..22}
do
  export gz=$(ls $prefix$i.*|grep -v log|wc -l)   
  export log=$(ls $prefix$i.*|grep .gz.log|wc -l)
  export ld=$(ls /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${i}_*_*.gz | wc -l)

#  if [ $gz -ne $log ] ;then 
  if [ $gz -ne $ld ] ;then 
    echo chr$i  gz=$gz log=$log ld=$ld
  fi
done
echo
}

## test in diff_MAX SNP per locus
for i in 1 3 5 7 10
do 
	cd ${path}/max_snp_$i
        pwd
       	echo $prefix/max_snp_${i}| tee -a $path/check_finemap.txt 
	check $prefix| tee -a $path/check_finemap.txt 
done


#path='/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_'
#path='/gpfs/commons/home/tlin/output/kunkle_all/finemap'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained'



  
