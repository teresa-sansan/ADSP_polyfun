#prefix='finemap_bellenguez_susie'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie'

prefix='finemap_bellenguez_all_2'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained'
agg='aggregrate.all.txt.gz'

echo $path

sig_level='1e-3'
for i in 1 3 5 7 10
do	
	## extract column name
	zcat ${path}/max_snp_$i/$agg | head -n 1 > ${path}/max_snp_$i/$prefix.extract_${sig_level}.csv 	
	echo start max_snp_$i ## extract those that reach GWAS significant level
	zcat ${path}/max_snp_$i/$agg | awk '$10< 1e-3 {print$0}' >> ${path}/max_snp_$i/$prefix.extract_${sig_level}.csv 
	#echo removing duplicated SNPs ....
	#cat ${path}/max_snp_$i/$prefix.extract_e-01.csv | awk '!a[$2]++' > ${path}/max_snp_$i/all_anno_extract_uniq_e-01.csv
	echo start zipping the file
	gzip ${path}/max_snp_$i/$prefix.extract_${sig_level}.csv 

done 


