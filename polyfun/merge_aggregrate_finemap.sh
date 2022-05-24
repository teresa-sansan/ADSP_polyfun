#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/'
path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/'

#summarystat=bellenguez
summary_stat='wightman'
for max_snp in 1 3 5 7 10
do 
  echo max_snp${max_snp}
  cd ${path}/max_snp_${max_snp}
  zcat chr11.aggregrate.all.txt.gz| head -n 1 > aggregate.all.txt
  for i in {1..22}
  do
     echo merging chr $i
     zcat chr${i}.aggregrate.all.txt.gz |grep -v CHR|cut -f 1-15 >> aggregate.all.txt 
  done
  echo
  echo zipping the file....
  gzip aggregate.all.txt 
  echo finished, total line = $( zcat aggregate.all.txt.gz | wc -l )
  

  zcat min_chr11.aggregrate.all.txt.gz| head -n 1 > agg_min_extract_1e-3.tsv
  echo 'extracting SNP with p < 0.001 ....'
  for chr in {1..22}
  do
  echo start aggregate chr $chr
### kunke & wightman p value are in 9th column, while bellengueez is in 10th
  zcat min_chr${chr}.aggregate.all.txt.gz| tail -n+2| awk '{if($9 <= 0.001) print$0}' >> agg_min_extract_1e-3.tsv
  done
  echo finished, total line = $( cat agg_min_extract_1e-3.tsv| wc -l )
done
