#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/'
path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
#path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824_no_partitions/'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'
#aggsummary_stat='bellenguez'
summary_stat='wightman'
#summary_stat='kunkle'
#summary_stat='jansen'

# bl no_ml all_anno all_enformer
#for anno in all_except_enformer bl no_ml_new
#for anno in susie all no_ml only_ml bl
for anno in bl
do
  for max_snp in 1  5 3 10
  do 
    echo max_snp${max_snp}
    cd ${path}/$anno/finemap/max_snp_${max_snp}
    #cd ${path}/finemap/max_snp_${max_snp}
    pwd
    zcat chr11.aggregate.all.txt.gz| head -n 1 > aggregate.all.txt
    #zcat chr11_adj_beta.aggregate.all.txt.gz| head -n 1 > adj_beta_aggregate.all.txt
    for i in {1..22}
    do
      echo merging chr $i
      zcat chr${i}.aggregate.all.txt.gz |tail -n+2  >> aggregate.all.txt
      #zcat chr${i}_adj_beta.aggregate.all.txt.gz |tail -n+2  >> adj_beta_aggregate.all.txt 

      #zcat chr${i}.aggregate.all.txt.gz |tail -n+2|awk '{if($9 <= 0.001) print$0}' >> aggregate.all.txt
      
      #echo zipping the file....
      #gzip aggregate.all.txt

      #echo finished, total line = $( zcat aggregate.all.txt.gz | wc -l )
    done

    ## see if want to create a smaller agg file.
      if true; then
      zcat chr11.aggregate.all.txt.gz| head -n 1 > agg_extract_1e-3.tsv
      #zcat chr11_adj_beta.aggregate.all.txt.gz| head -n 1 > adj_beta_agg_extract_1e-3.tsv
      echo "start extracting SNP with p < 1e-3 "
    ### kunke and wightman p value is in 9th column, while bellenguez and wightman are in 10th (sometime wightma is in 10)
      zcat chr${i}.aggregate.all.txt.gz|tail -n+2|awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
      cat aggregate.all.txt| tail -n+2| awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
      echo finished, total line = $( cat agg_extract_1e-3.tsv| wc -l )
    fi
    echo
done
done

