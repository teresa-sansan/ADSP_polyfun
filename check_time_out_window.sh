#cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap/max_snp_10
#touch missing_window.txt
for snp in 1 3 5 7 10
do
  cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_${snp}
  rm missing_window.txt
  touch missing_window.txt
  for chr in {1..22}
    do
    missing=0
     for i in finemap_bellenguez_all_2.${chr}.*.gz.log
       do
         export file=$(echo $i|awk -F '.' '{print $1 "." $2 "." $3 "." $4 ".gz"}')
         if  [ ! -f "$file" ]; then
            missing=$((++missing))
           if [ $missing -eq 1 ]; then
             echo | tee -a missing_window.txt 
             echo  chr $chr | tee -a missing_window.txt 
           fi
           echo $file  | tee -a missing_window.txt
         fi
       done
  
    done
done

