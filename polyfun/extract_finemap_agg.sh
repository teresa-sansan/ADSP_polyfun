prefix='agg_bellenguez'
path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap'
agg='aggregrate.all.txt.gz'

echo extracting finemapping result from $path ...

sig_level=1e-3
for i in 1 3 5 7 10
do	
	## extract column name
	rm ${path}/max_snp_$i/agg_bellenguez.extract_1e-3.tsv.gz
	echo original $(zcat ${path}/max_snp_$i/$agg | wc -l)  lines
	extract=${path}/max_snp_$i/$prefix.extract_${sig_level}.tsv
	echo $extract
	zcat ${path}/max_snp_$i/chr10.aggregrate.all.txt.gz | head -n 1 > $extract
	echo start max_snp_$i ## extract those that reach GWAS significant level
	zcat ${path}/max_snp_$i/$agg | awk '$9 < 1e-3 {print$0}' >> $extract
	echo extracted $(cat $extract | wc -l)  lines 
	echo start zipping the file
	gzip $extract 

done 


