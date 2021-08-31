name='/gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno'

zcat ${name}*.gz | awk '$9<1e-1 {print$0}'|uniq > ${name}_extract_e-01.csv 
echo start zipping the file
gzip ${name}_extract_e-01.csv 

zcat 

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
