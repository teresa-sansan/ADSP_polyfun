for i in finemap_max_snp_3_chr*.gz.log
do 
  gz=$(echo $i|cut -d '.' -f 1,2,3)
  if [ ! -f "$gz.gz" ]; then  
    if grep -Fq 'Estimating residual variance failed' $i ; then
       echo "negative residual variance"
       echo $gz
       echo ' ' 
    else
       echo "not finished", $gz
       echo ' '
    fi
  fi
done
