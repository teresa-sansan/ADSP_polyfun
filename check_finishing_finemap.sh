function check
{
for i in {1..22}
do
  export gz=$(ls finemap_bellenguez_all_2.$i.*|grep -v log|wc -l)   
  export log=$(ls finemap_bellenguez_all_2.$i.*|grep log|wc -l)
  if [ $gz -ne $log ] ;then 
    echo chr$i  gz=$gz log=$log
  fi
done
echo
}

path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap'
touch $path/check_finemap.txt

for i in 1 3 5 7 10
do 
	cd ${path}/max_snp_$i
       	echo bellenguez_all_2/max_snp_${i}| tee -a $path/check_finemap.txt 
	check| tee -a $path/check_finemap.txt 
done


#path='/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_'
#path='/gpfs/commons/home/tlin/output/kunkle_all/finemap'
#path='/gpfs/commons/home/tlin/output/kunkle_all_2/finemap'
