output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres'
sumstat='bellenguez'
for thres in 0.1 
do
    for anno in  susie
    do 
        # awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_pip${thres}_chr*.profile > $output/${anno}_pip${thres}.prs
        # echo wrote $output/${anno}_pip${thres}.prs
        awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_chr*.profile > $output/${anno}.prs
        echo wrote $output/${anno}.prs
    
    done
done

## baseline omics omics_dl 

#omics_dl_chr19.profile