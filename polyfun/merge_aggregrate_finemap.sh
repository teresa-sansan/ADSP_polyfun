#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/'


for max_snp in 1 3 5 7 10
do 
  echo max_snp${max_snp}
  cd ${path}/max_snp_${max_snp}
  zcat chr11.aggregrate.all.txt.gz| head -n 1 > aggregrate.all.txt
  for i in {1..22}
  do
     echo merging chr $i
     zcat chr${i}.aggregrate.all.txt.gz |grep -v CHR|cut -f 1-15 >> aggregrate.all.txt 
  done
  echo
  echo zipping the file....
  gzip aggregrate.all.txt 
  echo finished, total line = $( zcat aggregrate.all.txt.gz | wc -l )
done
