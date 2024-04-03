output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/plink_output/'
sumstat='bellenguez'
for thres in 0.1 
do 
    #awk '{ sum[$2]+=$4 } END {for (user in sum) print user, sum[user] }' $output/baseline_pip${thres}_chr*.profile > $output/variant_count_pip${thres}.txt
    #awk '{ sum[$2]+=$4 } END {for (user in sum) print user, sum[user] }' $output/baseline_chr*.profile > $output/variant_count_pip${thres}.txt
    
    awk '{ sum4[$2]+=$4; sum5[$2]+=$5} END { for (user in sum4) print user, sum4[user], sum5[user]}' $output/baseline_pip${thres}_chr*.profile > $output/variant_count_pip${thres}.txt
    awk '{ sum4[$2]+=$4; sum5[$2]+=$5} END { for (user in sum4) print user, sum4[user], sum5[user]}' $output/baseline_chr*.profile > $output/variant_count_no_thres.txt

    echo wrote $output/variant_count_pip${thres}.txt

done

## baseline omics omics_dl 

#omics_dl_chr19.profile