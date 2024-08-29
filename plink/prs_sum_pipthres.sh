output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/plink_output'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/allele_flip/'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/thres/plink_output/'
sumstat='bellenguez'
for thres in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
do
    for anno in susie baseline omics omics_dl
    do 
        # awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_chr*.profile > $output/${anno}.prs
        # echo wrote $output/${anno}.prs
        awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_chr*_pip_${thres}.profile > $output/${anno}_pip_${thres}.prs
        echo wrote  $output/${anno}_pip_${thres}.prs
    
    done
done

## baseline omics omics_dl  0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8