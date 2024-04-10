output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/p_thres'
sumstat='bellenguez'
for thres in 0.1 0.001 0.01 e-4 e-5 e-6
do
    for anno in susie
    do 
        # awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_pip${thres}_chr*.profile > $output/${anno}_pip${thres}.prs
        # echo wrote $output/${anno}_pip${thres}.prs
        awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_chr*.${thres}.profile > $output/${anno}_${thres}prs
        echo wrote  $output/${anno}_${thres}prs
    
    done
done

## baseline omics omics_dl 

#omics_dl_chr19.profile

# omics_chr5.e-4.profile  baseline omics_dl omics