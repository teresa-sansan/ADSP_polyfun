
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained'
prefix='finemap_bellenguez_all_2'
for i in 1 3 5 7 10
do	
	echo start max_snp_$i ## extract those that reach GWAS significant level
	zcat ${path}/max_snp_$i/$prefix*.gz | awk '$9<1e-5 {print$0}'|uniq > ${path}/max_snp_$i/$prefix.extract_e-01.csv 
	echo removing duplicated SNPs ....
	cat ${path}/max_snp_$i/$prefix.extract_e-01.csv | awk '!a[$2]++' > ${path}/max_snp_$i/all_anno_extract_uniq_e-01.csv
	echo start zipping the file
	#gzip ${path}/max_snp_$i/all_anno_extract_uniq_e-01.csv
 	break


done 


