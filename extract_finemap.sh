
path='/gpfs/commons/home/tlin/output/kunkle_all_2/finemap'
for i in 1 3 5 7 10
do	
	echo start max_snp_$i ## extract those that reach GWAS significant level
	zcat ${path}/max_snp_$i/all_anno*.gz | awk '$9<5*1e-8 {print$0}'|uniq > ${path}/max_snp_$i/anno_all.extract_e-01.csv 
	echo removing duplicated SNPs ....
	cat ${path}/max_snp_$i/anno_all.extract_e-01.csv | awk '!a[$2]++' > ${path}/max_snp_$i/all_anno_extract_uniq_e-01.csv
	echo start zipping the file
	gzip ${path}/max_snp_$i/all_anno_extract_uniq_e-01.csv
 

done 

#for i in H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9ac 
#do
#	cd ~/polyfun/output/bl_deepsea/$i/finemap
#	pwd
#	echo extracting $i
#	zcat finemap_bl_deepsea_${i}*.gz| awk '$9<1e-1 {print$0}'|uniq > finemap_bl_deepsea_${i}.extract_e-01.csv
#	echo start zipping the file
#	gzip finemap_bl_deepsea_${i}.extract_e-01.csv
#	
#done
