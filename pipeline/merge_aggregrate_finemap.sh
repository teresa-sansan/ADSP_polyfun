for max_snp in 1 3 5 7 
do 
  cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${max_snp}
  zcat chr1.aggregrate.all.txt.gz| head -n 1 > aggregrate.all.txt
  for i in {1..22}
  do
     echo merging chr $i
     zcat chr${i}.aggregrate.all.txt.gz |grep -v CHR >> aggregrate.all.txt
  done
  echo
  echo zipping the file....
  gzip aggregrate.all.txt 
  echo finished, total line = $( zcat aggregrate.all.txt | wc -l )
done