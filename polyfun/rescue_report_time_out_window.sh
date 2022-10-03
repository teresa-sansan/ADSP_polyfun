## Last edit: 10.03.22
## Purpose: Extract the regions that have gz.log but doesn't have gz.file due to time-out
## The reasons may be the following:
## (1) Those regions doesn't converge (2) Didn't request enough time.

#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
#prefix=finemap_bellenguez_all_2
path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/bl/finemap'
prefix='bl'
for snp in 1 5 10
do
  cd $path/max_snp_${snp}
  touch missing_window.txt
  for chr in {1..22}
    do
    missing=0
     for i in $prefix.chr${chr}.*.gz.log
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

