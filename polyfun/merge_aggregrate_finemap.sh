#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/new_susie/'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap/'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter/'
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/new_susie/'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/susie/'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/susie'
#path='/gpfs/commons/home/tlin/output/jansen/finemap/'
#path='/gpfs/commons/home/tlin/output/jansen/susie/'
#path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/all_anno/finemap/'
#path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/'

#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap/'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl/'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/'
path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/bl/finemap/'
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
     echo
     
     echo merging chr $i
     zcat chr${i}.aggregate.all.txt.gz |tail -n+2  >> aggregate.all.txt 
     #zcat chr${i}.aggregate.all.txt.gz |tail -n+2|awk '{if($9 <= 0.001) print$0}' >> aggregate.all.txt
     

    #echo zipping the file....
    #gzip aggregate.all.txt

    #echo finished, total line = $( zcat aggregate.all.txt.gz | wc -l )
    
    ## see if want to create a smaller agg file.
    if false; then
    zcat chr11.aggregate.all.txt.gz| head -n 1 > agg_extract_1e-3.tsv
    echo "start extracting SNP with p < 1e-3 "

  ### kunke  p value is in 9th column, while bellenguez and wightman are in 10th
    echo 'start extracting pvalue < 0.001'
    zcat chr${i}.aggregate.all.txt.gz|tail -n+2|awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
    #zcat aggregate.all.txt.gz| tail -n+2| awk '{if($10 <= 0.001) print$0}' >> agg_extract_1e-3.tsv
    echo finished, total line = $( cat agg_extract_1e-3.tsv| wc -l )
  fi
  done
done