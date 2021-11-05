#path='/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap'
#prefix='all_anno'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie'
#prefix=finemap_bellenguez_susie

path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
prefix=finemap_bellenguez_all_2
ld='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld'

for snp in 1 3 5 7 10
do
  cd $path/max_snp_${snp}
  if [ -f missing_window.txt ]; then rm missing_window.txt; fi
  touch missing_window.txt
  for chr in {1..22}
    do
    missing=0
     for i in $ld/chr${chr}_*_*.gz
       do  
         export file=$(echo $i| cut -d "_" -f 4,5 --output-delimiter='.' )
	 export file_name=$(echo $prefix.${chr}.${file})
         if  [ ! -f "$file_name" ]; then
            missing=$((++missing))
           if [ $missing -eq 1 ]; then
              echo | tee -a missing_window.txt 
              echo  chr $chr | tee -a missing_window.txt 
           fi
          echo $file_name  | tee -a missing_window.txt
         fi
       done
    done
done



#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/'
#prefix=finemap_bellenguez_all_2

