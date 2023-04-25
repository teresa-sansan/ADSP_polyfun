path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/no_ml_new/finemap/'
#summary_stat='bellenguez'
summary_stat='wightman'
#summary_stat='kunkle'
#summary_stat='jansen'

for max_snp in 1 5 10
do 
  echo max_snp${max_snp}
  cd ${path}/max_snp_${max_snp}
  pwd
  zcat chr11.aggregate.all.txt.gz| head -n 1 > aggregate.all.txt
  for i in {1..22}
  do
     echo merging chr $i
     zcat chr${i}.aggregate.all.txt.gz |tail -n+2  >> aggregate.all.txt 
     #zcat chr${i}.aggregate.all.txt.gz |tail -n+2|awk '{if($9 <= 0.001) print$0}' >> aggregate.all.txt
     
    #echo zipping the file....
    #gzip aggregate.all.txt

    #echo finished, total line = $( zcat aggregate.all.txt.gz | wc -l )
  done

  ## see if want to create a smaller agg file.
    if true; then
    zcat chr11.aggregate.all.txt.gz| head -n 1 > agg_extract_1e-3.tsv
    echo "start extracting SNP with p < 1e-3 "
  ### kunke  p value is in 9th column, while bellenguez and wightman are in 10th
    #zcat chr${i}.aggregate.all.txt.gz|tail -n+2|awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
    cat aggregate.all.txt| tail -n+2| awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
    echo finished, total line = $( cat agg_extract_1e-3.tsv| wc -l )
  fi
  echo
done

