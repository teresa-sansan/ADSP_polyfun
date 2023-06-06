path='/gpfs/commons/home/tlin/output/bellenguez/new_anno/'
prefix='bellenguez.'

## create a function to check if its converged
function check
{
for i in {1..22}
do
  export gz=$(ls $prefix$i.*|grep -v log|wc -l)   
  export log=$(ls $prefix$i.*|grep .gz.log|wc -l)
  export ld=$(ls /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${i}_*_*.gz | wc -l)


  if [ $gz -ne $ld ] ;then 
    echo chr$i  gz=$gz log=$log ld=$ld
  fi
  #if [ $gz -eq $ld ] ;then
  # echo checked chr$i
  #fi
done
echo
}

## test in diff_MAX SNP per locus
for anno in no_ml bl
do
if [ -f $path/$anno/check_finemap.txt ]; then rm $path/$anno/check_finemap.txt; fi #remove pre-exist file
touch $path/$anno/check_finemap.txt   ##create file

 for i in 2 3
 do 
	cd ${path}/$anno/finemap/max_snp_$i
        pwd
       	echo $prefix/max_snp_${i}| tee -a $path/$anno/check_finemap.txt 
	check $prefix| tee -a $path/$anno/check_finemap.txt 
 done
done





  
