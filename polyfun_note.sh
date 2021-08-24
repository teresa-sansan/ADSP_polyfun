##Following the tutorial on polyfun wiki
#Apporch1 -- precomputed prior

##testing data
##pwd ##~/polyfun

mkdir -p output 
conda activate polyfun
python extract_snpvar.py --sumstats example_data/snps_to_finemap.txt.gz --out output/snps_with_│
var.gz
cd output/
zcat snps_with_var.gz| head


##AD data
##pwd ##~/polyfun

mkdir try_Kunkle
python extract_snpvar.py --sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary│
_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet --out try_Kunkle/kunkle_with_var.gz 

## There are missing SNPs in the file, might miss some causal SNPs
#outputfile: ~/polyfun/test_kunkle/kunckle_with_var.gz.miss.gz

zcat kunkle_with_var.gz.miss.gz| wc   #1340939 9385299 85944550

## Add "allow-missing" flag. Notice that we may lose causal SNP
python extract_snpvar.py --sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stat/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet --out try_Kunkle/kunkle_allowmissing.gz --allow-missing


##check kunckle.parquet file
module load parquet-tools

parquet-tools inspect /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheim│ers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet | head -n 30 >  try_Kunkle/parquethead30.txt   
## Appearing broken pipe error. 


#Apporch2 -- computing prior using L2 extension of S-LDSC
## first I need to install each pack individually. So I removed the environment and reinstall the polyfun environment. 
## toy example

mkdir -p output

python polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix output/testrun \
    --sumstats example_data/sumstats.parquet \
    --ref-ld-chr example_data/annotations. \
    --w-ld-chr example_data/weights.


##



extract file (CHR SNP SP)
for f in baselineLD.*.annot.gz; do zcat $f| cut -f 1,2,3 |gzip> ~/polyfun/data_backup/anno.$f ; done 


for f in brain_atac_seq_chr*l2.ldscore.gz; do cp $f ~/polyfun/data_backup/anno.$f; done

for f in brain_atac_seq_chr*.annot.gz ; do  zcat $f | cut -f 1,3,2,5,6,8 |gzip > ~/polyfun/data_backup/anno.$f ; done    


python  polyfun.py --compute-h2-L2 --output-prefix output/testrun --sumstats data/kunkle_2019/AD_Kunkle_etal_Stage1.parquet --ref-ld-chr data_backup/brain_atac_polyfun/brain_atac_seq_chr  --w-ld-chr /gpfs/commons/groups/knowles_lab/data/ldsc/weights_hm3_no_hla/weights.  

for i in brain_atac_seq_chr*.annot.gz        
  do                                                                   
  zcat $i| cut -f 1,2,3,5,6,8| gzip > ~/polyfun/data/$i                               
  done                                                 


##success
for f in brain_atac_seq_chr*.annot.gz; do;zcat $f|cut -f 1,2,3,5,6,8| gzip > ~/polyfun/data/$f;done

#edit ldscore file #changed microglia_atacL2 to microglia_atac
for f in brain_atac_seq_chr*.l2.ldscore.gz
do
zcat $f | cut -f 1,2,3,4,5,7,| sed 's.microglia_atacL2/microglia_atac'| gzip -c > ~/polyfun/fata/$f
done
