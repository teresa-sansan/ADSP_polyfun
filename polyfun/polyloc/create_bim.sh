path='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/UKB_snp_list'
cd $path
for i in {1..22}
do
    echo start chr $i ...
    echo CHR$'\t'SNP$'\t'CM$'\t'BP$'\t'A1$'\t'A2 > UKB_SNP_list_chr${i}.bim
    zcat UKB_SNP_list_chr${i}.gz|awk '{print $1,$2,'0',$3,$4,$5}'|tr ' ' '\t'| tail -n+2 >> UKB_SNP_list_chr${i}.bim
done