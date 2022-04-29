echo max_snp = 3
echo 'max_snp=3' > not_converge.txt
for i in {1..3}
do 
 echo chr$i
 total=$(ls|grep max_snp_3_chr$i| grep log|wc -l)
 echo total = $total
 not_converge=$(cat finemap_max_snp_3_chr$i.*.gz.log| grep IBSS| wc -l)
 negative=$(cat finemap_max_snp_3_chr$i.*.gz.log| grep 'the estimated value is negative'| wc -l) 
 echo not converge = $not_converge , negative residual variance = $negative
 if [ $not_converge != 0 ]; then
 for log in finemap_max_snp_3_chr${i}.*.gz.log
      do if grep -Fq 'IBSS' $log; then  echo $log >> not_converge.txt; fi; done
 fi


done

if true;then
for max_snp in 5 7 10
do
echo
echo max_snp = $max_snp 
 for i in 1 2; 
  do echo chr$i
  total=$(ls|grep max_snp_${max_snp}_chr${i}| grep log|wc -l)
  echo total = $total 
  echo not converge = $(cat finemap_max_snp_${max_snp}_chr$i.*.gz.log| grep IBSS| wc -l) , negative residual variance = $(cat finemap_max_snp_${max_snp}_chr$i.*.gz.log| grep negative| wc -l)

  if grep -Fq IBSS finemap_max_snp_${max_snp}_chr${i}.*.gz.log; then
     echo >> not_converge.txt
     echo max_snp_${max_snp} >> not_converge.txt
     echo IBSS >> not_converge.txt
     echo finemap_max_snp_${max_snp}_chr${i}.*.gz.log >> not_converge.txt;
  fi

#  if grep -Fq 'the estimated value is negative' finemap_max_snp_${max_snp}_chr${i}.*.gz.log; then
#    echo ' ' >> not_converge.txt
#    echo negative_residual_variance >> not_converge.txt
#    echo finemap_max_snp_${max_snp}_chr${i}.*.gz.log >> not_converge.txt
#  fi

  done
done
fi
cat not_converge.txt | sed s/'finemap_'//g | sed s/'.gz.log'//g| sed s/' '/'\n'/g > not_converge2.txt 
mv not_converge2.txt not_converge.txt
#cat not_converge.txt
