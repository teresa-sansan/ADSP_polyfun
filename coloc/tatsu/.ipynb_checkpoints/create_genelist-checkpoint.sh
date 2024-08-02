cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia

if [ -e genelist ];then
    rm genelist
else
    touch genelist
fi

for chr in {1..22}
do
echo chr $chr
    cat microglia_eqtl_chr${chr}.tsv | awk '$31 < 1e-05 {print $1}'|uniq |sort|uniq > chr${chr}_gene 
    cat chr${chr}_gene | wc 
    cat chr${chr}_gene >> genelist
done