#cd /gpfs/commons/home/tlin/output/others/qc
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/plink_biallelic/qc_maf_0.1

for i in {1..22}
do


file=target_QC_chr${i}.log
writefile=chr${i}.info

echo chr${i} > $writefile
echo ' '  >> $writefile 
echo 'before QC:' >>$writefile
cat $file| grep "variants loaded from .bim file."| head -1 >> $writefile
cat $file|grep "loaded from .fam."| head -1 >> $writefile
echo ' '  >> $writefile 
cat $file| grep "removed"  >> $writefile
echo ' ' >> $writefile  
cat $file| grep 'Total genotyping rate'  >> $writefile 
echo ' '  >> $writefile 
cat $file| grep "pass filter"  >> $writefile 
echo ' '  >> $writefile 


done

