cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/
#touch gene_pos.tsv
#touch gene_range.tsv
rm gene_tss.tsv
touch gene_tss.tsv
cat genelist | xargs -P 5 -I {} bash -c '
    gene="{}"
    info=$(python /gpfs/commons/home/tlin/script/coloc/get_parquet_tss.py "$gene")
    read chr tss_pos <<< "$info"
    echo "$gene $chr $tss_pos" | tee -a gene_tss.tsv
'
