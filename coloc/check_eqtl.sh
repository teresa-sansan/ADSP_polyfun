#rm /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene3.txt
#touch /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene3.txt
rm /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss_rerun_gene2.txt
touch /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss_rerun_gene2.txt
cat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/genelist | while IFS=$'\t' read -r gene
#tail -n +40 /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene.txt | while IFS=$'\t' read -r gene
do
    #if [ ! -e /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_z2_chr__${gene}.rds ]; then
    if [ ! -e /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_chr*_${gene}.rds ]; then
        #echo running $gene
        ##bash eqtl_finemap.sh $gene
        #echo $gene >> /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene3.txt      
        echo $gene >> /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss_rerun_gene2.txt
    fi      
done


# echo start running eqtl
# #cat /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene2.txt |
# #tail -n +500 /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/rerun_gene3.txt | \
# cat /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss_rerun_gene.txt
# xargs -P 5 -I {} bash -c '
#     gene="{}"
# #   if [ ! -e /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl/eqtl_z2_chr__${gene}.rds ]; then
#     if [ ! -e /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss_oct/eqtl_tss/eqtl_chr*_${gene}.rds ]; then
#         echo running $gene
#         bash eqtl_finemap.sh $gene
#     else
#         echo $gene is finished
#     fi
# '